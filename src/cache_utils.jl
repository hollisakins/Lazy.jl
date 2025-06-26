"""
Template Grid Caching Utilities for Lazy.jl

Provides efficient caching of template grids to speed up subsequent runs.
Cache is stored in HDF5 format with automatic cleanup and invalidation.
"""

using HDF5
using Glob

const CACHE_DIR = joinpath(@__DIR__, "cache")
const MAX_CACHE_SIZE_GB = 5.0
const MAX_CACHE_SIZE_BYTES = Int64(MAX_CACHE_SIZE_GB * 1024^3)

"""
    ensure_cache_dir()

Create cache directory if it doesn't exist.
"""
function ensure_cache_dir()
    if !isdir(CACHE_DIR)
        mkpath(CACHE_DIR)
    end
end

"""
    generate_cache_key(templates::Vector{String}, zgrid::Vector{Float64}, 
                      bands::Vector{String}, igm_model::String, 
                      template_error::String, template_error_scale::Float64)::String

Generate a unique cache key based on all parameters that affect template grid generation.
Returns an 8-character hex string.
"""
function generate_cache_key(templates::Vector{String}, zgrid::Vector{Float64}, 
                           bands::Vector{String}, igm_model::String, 
                           template_error::String, template_error_scale::Float64)::String
    
    # Hash template files and their modification times
    template_mtimes = [stat(t).mtime for t in templates]
    template_hash = hash((templates, template_mtimes))
    
    # Hash redshift grid parameters
    grid_hash = hash((zgrid[1], zgrid[end], length(zgrid)))
    
    # Hash bands (order matters for template grid indexing)
    bands_hash = hash(bands)
    
    # Hash IGM model and its modification time
    igm_file = joinpath(@__DIR__, "igm_data", igm_model * ".hdf5")
    igm_mtime = isfile(igm_file) ? stat(igm_file).mtime : 0.0
    igm_hash = hash((igm_model, igm_mtime))
    
    # Hash template error parameters
    template_error_file = joinpath(@__DIR__, "templates", "template_error", template_error * ".dat")
    template_error_mtime = isfile(template_error_file) ? stat(template_error_file).mtime : 0.0
    error_hash = hash((template_error, template_error_scale, template_error_mtime))
    
    # Combine all hashes
    combined_hash = hash((template_hash, grid_hash, bands_hash, igm_hash, error_hash))
    
    return string(combined_hash, base=16)[1:8]
end

"""
    get_cache_filename(cache_key::String)::String

Get the full path to a cache file for the given cache key.
"""
function get_cache_filename(cache_key::String)::String
    return joinpath(CACHE_DIR, "templgrid_$(cache_key).h5")
end

"""
    cache_exists(cache_key::String)::Bool

Check if a valid cache file exists for the given cache key.
"""
function cache_exists(cache_key::String)::Bool
    cache_file = get_cache_filename(cache_key)
    return isfile(cache_file)
end

"""
    save_template_cache(templgrid::Array{Float64,3}, template_error_grid::Matrix{Float64},
                       metadata::Dict, cache_key::String)

Save template grid and associated data to cache file.
"""
function save_template_cache(templgrid::Array{Float64,3}, template_error_grid::Matrix{Float64},
                           metadata::Dict, cache_key::String)
    ensure_cache_dir()
    cache_file = get_cache_filename(cache_key)
    
    try
        h5open(cache_file, "w") do file
            # Save main template grid
            file["templgrid"] = templgrid
            
            # Save template error grid
            file["template_error_grid"] = template_error_grid
            
            # Save metadata
            grp = create_group(file, "metadata")
            for (key, value) in metadata
                if value isa Vector{String}
                    grp[key] = value
                elseif value isa String || value isa Number
                    grp[key] = value
                elseif value isa Vector{<:Number}
                    grp[key] = value
                else
                    # Convert to string for complex types
                    grp[key] = string(value)
                end
            end
            
            # Add timestamp
            grp["created"] = time()
            grp["cache_key"] = cache_key
        end
        
        # Update access time for LRU cleanup
        try
            touch(cache_file)
        catch
            # Fallback: just ensure file exists (touch may fail on some systems)
        end
        
        # Cleanup cache if needed
        cleanup_cache()
        
    catch e
        @warn "💾 Failed to save template cache: $e"
        # Remove partially written file
        if isfile(cache_file)
            rm(cache_file, force=true)
        end
    end
end

"""
    load_template_cache(cache_key::String)::Union{Tuple{Array{Float64,3}, Matrix{Float64}}, Nothing}

Load template grid from cache. Returns (templgrid, template_error_grid) if successful, nothing if failed.
"""
function load_template_cache(cache_key::String)::Union{Tuple{Array{Float64,3}, Matrix{Float64}}, Nothing}
    cache_file = get_cache_filename(cache_key)
    
    if !isfile(cache_file)
        return nothing
    end
    
    try
        templgrid = nothing
        template_error_grid = nothing
        
        h5open(cache_file, "r") do file
            templgrid = read(file, "templgrid")
            template_error_grid = read(file, "template_error_grid")
        end
        
        # Update access time for LRU
        try
            touch(cache_file)
        catch
            # Fallback: touch may fail on some systems
        end
        
        return (templgrid, template_error_grid)
        
    catch e
        @warn "💾 Failed to load template cache $(cache_key): $e"
        # Remove corrupted cache file
        rm(cache_file, force=true)
        return nothing
    end
end

"""
    get_cache_info()::Dict

Get information about the current cache state.
"""
function get_cache_info()::Dict
    ensure_cache_dir()
    cache_files = glob("templgrid_*.h5", CACHE_DIR)
    
    total_size = isempty(cache_files) ? 0 : sum(stat(f).size for f in cache_files)
    file_count = length(cache_files)
    
    return Dict(
        "file_count" => file_count,
        "total_size_bytes" => total_size,
        "total_size_gb" => total_size / (1024^3),
        "max_size_gb" => MAX_CACHE_SIZE_GB,
        "directory" => CACHE_DIR
    )
end

"""
    cleanup_cache()

Remove old cache files to maintain size limit using LRU (least recently used) policy.
"""
function cleanup_cache()
    ensure_cache_dir()
    cache_files = glob("templgrid_*.h5", CACHE_DIR)
    
    if isempty(cache_files)
        return
    end
    
    # Calculate total size with fallback for atime access issues
    file_info = []
    for f in cache_files
        try
            file_stat = stat(f)
            file_size = file_stat.size
            
            # Try to access atime, fallback to mtime if not available (macOS/Julia compatibility)
            access_time = try
                file_stat.atime
            catch
                file_stat.mtime
            end
            
            push!(file_info, (f, file_size, access_time))
        catch e
            @warn "🧹 Could not process cache file $f: $e"
            # Skip this file if we can't stat it
        end
    end
    
    if isempty(file_info)
        return  # No valid cache files to process
    end
    
    total_size = sum(info[2] for info in file_info)
    
    if total_size <= MAX_CACHE_SIZE_BYTES
        return  # No cleanup needed
    end
    
    # Sort by access time (oldest first)
    sort!(file_info, by=x -> x[3])
    
    # Remove files until under limit
    removed_count = 0
    for (file_path, file_size, _) in file_info
        if total_size <= MAX_CACHE_SIZE_BYTES
            break
        end
        
        try
            rm(file_path)
            total_size -= file_size
            removed_count += 1
        catch e
            @warn "🧹 Failed to remove cache file $(file_path): $e"
        end
    end
    
    if removed_count > 0
        println("🧹 Cache cleanup: removed $removed_count old cache files")
    end
end

"""
    clear_cache()

Remove all cache files.
"""
function clear_cache()
    ensure_cache_dir()
    cache_files = glob("templgrid_*.h5", CACHE_DIR)
    
    removed_count = 0
    freed_bytes = 0
    
    for cache_file in cache_files
        try
            file_size = stat(cache_file).size
            rm(cache_file)
            removed_count += 1
            freed_bytes += file_size
        catch e
            @warn "🧹 Failed to remove cache file $(cache_file): $e"
        end
    end
    
    freed_gb = freed_bytes / (1024^3)
    return (removed_count, freed_gb)
end