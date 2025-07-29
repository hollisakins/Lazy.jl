"""
HDF5 Streaming functionality for Lazy.jl

Provides streaming write capabilities, resume functionality, and HDF5-FITS conversion.
"""

using HDF5
using OrderedCollections
using FITSIO
using TOML

"""
    store_param_metadata(parent_group::HDF5.Group, param::Dict)

Store parameter file contents as structured HDF5 groups and datasets.
This avoids issues with TOML.print() returning Nothing.
"""
function store_param_metadata(parent_group::HDF5.Group, param::Dict)
    param_group = create_group(parent_group, "parameters")
    
    for (section_name, section_data) in param
        if isa(section_data, Dict)
            # Create group for each section (io, fitting, runtime, translate)
            section_group = create_group(param_group, section_name)
            
            for (key, value) in section_data
                # Store individual parameters as datasets
                try
                    if isa(value, Dict)
                        # Handle nested dictionaries (like in translate section)
                        nested_group = create_group(section_group, key)
                        for (nested_key, nested_value) in value
                            nested_group[nested_key] = string(nested_value)
                        end
                    elseif isa(value, Vector)
                        # Handle arrays (like template file lists)
                        section_group[key] = string.(value)
                    else
                        # Handle scalar values
                        section_group[key] = string(value)
                    end
                catch e
                    # If storage fails, store as string representation
                    println("Warning: Could not store parameter $section_name.$key as structured data, storing as string")
                    section_group[key] = string(value)
                end
            end
        else
            # Handle top-level non-dict values
            param_group[section_name] = string(section_data)
        end
    end
end

"""
    create_hdf5_work_file(filename::String, param::Dict, nobj::Int, nz::Int, 
                          nband::Int, ntempl::Int, bands::Vector{String}, 
                          zgrid::Vector{Float64}, templates::Vector{String})

Create HDF5 work file with proper structure and metadata for streaming writes.
"""
function create_hdf5_work_file(filename::String, param::Dict, nobj::Int, nz::Int, 
                              nband::Int, ntempl::Int, bands::Vector{String}, 
                              zgrid::Vector{Float64}, templates::Vector{String})
    
    h5open(filename, "w") do file
        # Create metadata group
        g_meta = create_group(file, "metadata")
        
        # Store parameter file as structured HDF5 groups instead of TOML string
        store_param_metadata(g_meta, param)
        
        # Store version safely
        version_str = try
            version
        catch
            "unknown"
        end
        g_meta["version"] = version_str
        
        g_meta["start_time"] = time()
        g_meta["nobj"] = nobj
        g_meta["nz"] = nz
        g_meta["nband"] = nband
        g_meta["ntempl"] = ntempl
        g_meta["bands"] = bands
        g_meta["zgrid"] = zgrid
        g_meta["templates"] = templates
        g_meta["last_processed_object"] = 0
        g_meta["last_update_time"] = time()
        g_meta["complete"] = false
        
        # Create results group with pre-allocated datasets
        g_results = create_group(file, "results")
        create_dataset(g_results, "ID", Int32, (nobj,))
        create_dataset(g_results, "zbest", Float64, (nobj,))
        create_dataset(g_results, "chi2best", Float64, (nobj,))
        create_dataset(g_results, "z_l95", Float64, (nobj,))
        create_dataset(g_results, "z_l68", Float64, (nobj,))
        create_dataset(g_results, "z_med", Float64, (nobj,))
        create_dataset(g_results, "z_u68", Float64, (nobj,))
        create_dataset(g_results, "z_u95", Float64, (nobj,))
        create_dataset(g_results, "coeffs", Float64, (nobj, ntempl))
        
        # Create photometry group for best-fit model
        g_phot = create_group(file, "photometry")
        for (i, band) in enumerate(bands)
            create_dataset(g_phot, band, Float64, (nobj,))
        end
        
        # Optional P(z) group with compression
        output_pz = get(param["io"], "output_pz", false)
        if output_pz
            g_pz = create_group(file, "pz")
            # Use chunking and compression for large chi2 grid
            chunk_size = min(1000, nobj)
            create_dataset(g_pz, "chi2grid", Float32, (nobj, nz), 
                          chunk=(chunk_size, nz), compress=3)
        end
        
        # Optional templates group (datasets created when template data is written)
        output_templ = get(param["io"], "output_templates", false)
        if output_templ
            g_templates = create_group(file, "templates")
            # Template datasets will be created dynamically when template data is written
            # This avoids pre-allocating without knowing wavelength grid size
        end
    end
end

"""
    check_resume_file(filename::String)

Check if a work file exists and is resumable. Returns (can_resume, last_object).
"""
function check_resume_file(filename::String)
    if !isfile(filename)
        return (false, 0)
    end
    
    try
        h5open(filename, "r") do file
            if !haskey(file, "metadata")
                return (false, 0)
            end
            
            meta = file["metadata"]
            complete = read(meta, "complete")
            if complete
                return (false, 0)  # Already complete
            end
            
            last_obj = read(meta, "last_processed_object")
            return (true, last_obj)
        end
    catch
        return (false, 0)  # Corrupted file
    end
end

"""
    write_chunk_results(filename::String, chunk_start::Int, chunk_end::Int,
                       IDs::Vector{<:Integer}, zbest::Vector{Float64}, 
                       chi2best::Vector{Float64}, coeffsbest::Matrix{Float64},
                       z_l95::Vector{Float64}, z_l68::Vector{Float64},
                       z_med::Vector{Float64}, z_u68::Vector{Float64},
                       z_u95::Vector{Float64}, photobest::Matrix{Float64},
                       bands::Vector{String}, chi2grid::Union{Nothing,Matrix{Float32}}=nothing)

Write a chunk of results to the HDF5 work file.
"""
function write_chunk_results(filename::String, chunk_start::Int, chunk_end::Int,
                           IDs::Vector{<:Integer}, zbest::Vector{Float64}, 
                           chi2best::Vector{Float64}, coeffsbest::Matrix{Float64},
                           z_l95::Vector{Float64}, z_l68::Vector{Float64},
                           z_med::Vector{Float64}, z_u68::Vector{Float64},
                           z_u95::Vector{Float64}, photobest::Matrix{Float64},
                           bands::Vector{String}, chi2grid::Union{Nothing,Matrix{Float32}}=nothing)
    
    try
        h5open(filename, "r+") do file
            # Write results (convert IDs to Int32 to match dataset type)
            file["results/ID"][chunk_start:chunk_end] = Int32.(IDs)
            file["results/zbest"][chunk_start:chunk_end] = zbest
            file["results/chi2best"][chunk_start:chunk_end] = chi2best
            file["results/z_l95"][chunk_start:chunk_end] = z_l95
            file["results/z_l68"][chunk_start:chunk_end] = z_l68
            file["results/z_med"][chunk_start:chunk_end] = z_med
            file["results/z_u68"][chunk_start:chunk_end] = z_u68
            file["results/z_u95"][chunk_start:chunk_end] = z_u95
            file["results/coeffs"][chunk_start:chunk_end, :] = coeffsbest
            
            # Write photometry
            for (i, band) in enumerate(bands)
                file["photometry/$band"][chunk_start:chunk_end] = photobest[:, i]
            end
            
            # Write P(z) if provided
            if chi2grid !== nothing && haskey(file, "pz/chi2grid")
                file["pz/chi2grid"][chunk_start:chunk_end, :] = chi2grid
            end
            
            # Update metadata (handle existing vs new datasets)
            meta_group = file["metadata"]
            
            # Update existing dataset
            if haskey(meta_group, "last_processed_object")
                delete_object(meta_group, "last_processed_object")
            end
            meta_group["last_processed_object"] = chunk_end
            
            # Create or update last_update_time
            if haskey(meta_group, "last_update_time")
                delete_object(meta_group, "last_update_time")
            end
            meta_group["last_update_time"] = time()
            
            # Force flush to disk
            flush(file)
        end
    catch e
        # Provide safe error message without printing large arrays
        println("❌ Error writing chunk results to $filename:")
        println("   Chunk: objects $chunk_start to $chunk_end")
        println("   Array sizes: IDs=$(length(IDs)), zbest=$(length(zbest))")
        
        # Safe error message handling
        error_msg = if hasfield(typeof(e), :msg)
            e.msg
        else
            string(e)
        end
        
        println("   Error: $(typeof(e)) - $error_msg")
        rethrow(e)
    end
end

"""
    write_template_data(work_file::String, templgrid_spectral, wavelength_grid, zgrid, templates)

Write template spectral data to HDF5 work file.
templgrid_spectral has shape (ntempl, nz, nwav) and is converted to (nwav, nz) per template for storage.
"""
function write_template_data(work_file::String, templgrid_spectral, wavelength_grid, zgrid, templates)
    h5open(work_file, "r+") do file
        g_templates = file["templates"]
        
        # Store wavelength grid (common to all templates)
        g_templates["wavelength"] = wavelength_grid
        
        # Store each template's spectral data
        # templgrid_spectral[i,:,:] has shape (nz, nwav)
        # We need to transpose to (nwav, nz) to match original format
        for (i, template) in enumerate(templates)
            template_name = basename(template)
            template_name = splitext(template_name)[1]  # Remove extension regardless of type
            
            # Convert from (nz, nwav) to (nwav, nz) to match original format
            template_data = transpose(templgrid_spectral[i, :, :])
            # Convert transpose to regular array for HDF5 compatibility
            g_templates[template_name] = collect(template_data)
        end
    end
end

"""
    finalize_hdf5_work_file(filename::String)

Mark HDF5 work file as complete.
"""
function finalize_hdf5_work_file(filename::String)
    h5open(filename, "r+") do file
        meta_group = file["metadata"]
        
        # Update complete status
        if haskey(meta_group, "complete")
            delete_object(meta_group, "complete")
        end
        meta_group["complete"] = true
        
        # Update end time
        if haskey(meta_group, "end_time")
            delete_object(meta_group, "end_time")
        end
        meta_group["end_time"] = time()
    end
end

"""
    convert_hdf5_to_fits(hdf5_file::String, fits_file::String)

Convert completed HDF5 work file to FITS format.
"""
function convert_hdf5_to_fits(hdf5_file::String, fits_file::String)
    h5open(hdf5_file, "r") do h5f
        # Read metadata
        meta = h5f["metadata"]
        bands = read(meta, "bands")
        nobj = read(meta, "nobj")
        
        # Read results
        IDs = read(h5f["results/ID"])
        zbest = read(h5f["results/zbest"])
        chi2best = read(h5f["results/chi2best"])
        z_l95 = read(h5f["results/z_l95"])
        z_l68 = read(h5f["results/z_l68"])
        z_med = read(h5f["results/z_med"])
        z_u68 = read(h5f["results/z_u68"])
        z_u95 = read(h5f["results/z_u95"])
        coeffsbest = read(h5f["results/coeffs"])
        
        # Prepare data for FITS
        data = OrderedDict{String, Dict{String, Any}}()
        data["ID"] = Dict("format" => "K", "data" => IDs)
        data["z_best"] = Dict("format" => "E", "data" => zbest)
        data["chi2"] = Dict("format" => "E", "data" => chi2best)
        data["z_l95"] = Dict("format" => "E", "data" => z_l95)
        data["z_l68"] = Dict("format" => "E", "data" => z_l68)
        data["z_med"] = Dict("format" => "E", "data" => z_med)
        data["z_u68"] = Dict("format" => "E", "data" => z_u68)
        data["z_u95"] = Dict("format" => "E", "data" => z_u95)
        
        # Add photometry
        for (i, band) in enumerate(bands)
            photobest_band = read(h5f["photometry/$band"])
            data[band] = Dict("format" => "E", "unit" => "fnu", "data" => photobest_band)
        end
        
        ntempl = size(coeffsbest, 2)
        data["coeffs"] = Dict("format" => "$(ntempl)E", "data" => coeffsbest)
        
        # Write SUMMARY extension
        write_data(fits_file, data, "SUMMARY")
        
        # Write P(z) if it exists
        if haskey(h5f, "pz/chi2grid")
            zgrid = read(meta, "zgrid")
            chi2grid = read(h5f["pz/chi2grid"])
            
            # Convert chi2 to P(z)
            pz = exp.(-0.5 * chi2grid)
            pz = pz ./ trapz(zgrid, pz)
            
            # Format for FITS
            temp_pz = vcat(transpose(zgrid), pz)
            temp_IDs = vcat([-1], IDs)
            
            nz = length(zgrid)
            data_pz = OrderedDict{String, Dict{String, Any}}()
            data_pz["ID"] = Dict("format" => "K", "data" => temp_IDs)
            data_pz["Pz"] = Dict("format" => "$(nz)E", "data" => temp_pz)
            
            write_data(fits_file, data_pz, "PZ")
        end
        
        # Write TEMPL extension if it exists
        if haskey(h5f, "templates")
            zgrid = read(meta, "zgrid")
            templates = read(meta, "templates")
            wavelength_grid = read(h5f["templates/wavelength"])
            
            nwav = length(wavelength_grid)
            
            # Prepare data for FITS (matching original format exactly)
            data_templ = OrderedDict{String, Dict{String, Any}}()
            temp_zgrid = vcat([-1], zgrid)
            data_templ["z"] = Dict("format" => "E", "data" => temp_zgrid)
            
            for template in templates
                template_name = basename(template)
                template_name = splitext(template_name)[1]
                
                # Read template data (stored as (nwav, nz) format)
                template_data = read(h5f["templates/$template_name"])
                
                # Format for FITS output (matching original exactly)
                # Need to transpose template_data back to (nz, nwav) for FITS format
                # FITS format: (nz+1, nwav) where first row is wavelength, then nz rows of template flux
                template_data_transposed = transpose(template_data)  # Convert (nwav, nz) -> (nz, nwav)
                temp_template_data = vcat(transpose(wavelength_grid), template_data_transposed)
                data_templ[template_name] = Dict("format" => "$(nwav)E", "data" => temp_template_data)
            end
            
            write_data(fits_file, data_templ, "TEMPL")
        end
    end
end

"""
    get_optimal_chunk_size(nobj::Int, nz::Int)

Calculate optimal chunk size based on default memory target (0.5 GB).
Maintained for backwards compatibility.
"""
function get_optimal_chunk_size(nobj::Int, nz::Int)
    return get_chunk_size_for_memory_target(nobj, nz, 0.5)
end

"""
    get_chunk_size_for_memory_target(nobj::Int, nz::Int, target_memory_gb::Float64)

Calculate chunk size to achieve target memory usage per chunk.
"""
function get_chunk_size_for_memory_target(nobj::Int, nz::Int, target_memory_gb::Float64)
    bytes_per_object = nz * 4  # Float32 chi2 values
    target_bytes = target_memory_gb * 1024^3
    objects_per_chunk = floor(Int, target_bytes / bytes_per_object)
    
    # Clamp to reasonable range (min 1000 objects, max total objects)
    return clamp(objects_per_chunk, 1_000, nobj)
end

"""
    prompt_complete_run(work_file::String, last_obj::Int, nobj::Int)

Prompt user when a run is already complete.
Returns: action choice ("complete", "overwrite", "keep")
"""
function prompt_complete_run(work_file::String, last_obj::Int, nobj::Int)::String
    # Get file modification time
    mtime = stat(work_file).mtime
    time_ago = time() - mtime
    time_str = if time_ago < 3600
        "$(round(Int, time_ago / 60)) minutes ago"
    elseif time_ago < 86400
        "$(round(time_ago / 3600, digits=1)) hours ago"
    else
        "$(round(Int, time_ago / 86400)) days ago"
    end
    
    println("🧠 Found complete run (100.0% complete, $(format_number(last_obj))/$(format_number(nobj)) objects)")
    println("   Last updated: $time_str")
    println("   Options:")
    println("   [C]omplete - Mark as complete and create final output")
    println("   [O]verwrite - Delete work file and start fresh")
    println("   [K]eep - Keep work file and exit")
    print("   Choice [C/o/k]: ")
    
    response = lowercase(strip(readline()))
    
    if response == "o"
        return "overwrite"
    elseif response == "k"
        return "keep"
    else  # Default to complete
        return "complete"
    end
end

"""
    prompt_resume(work_file::String, last_obj::Int, nobj::Int)

Prompt user whether to resume from checkpoint.
Returns: (action, should_resume) where action is "resume", "overwrite", "complete", or "keep"
"""
function prompt_resume(work_file::String, last_obj::Int, nobj::Int)::Tuple{String, Bool}
    # Check if the run is actually complete
    if last_obj >= nobj
        action = prompt_complete_run(work_file, last_obj, nobj)
        return (action, false)  # No resuming needed for complete runs
    end
    
    percent_complete = round(100 * last_obj / nobj, digits=1)
    
    # Get file modification time
    mtime = stat(work_file).mtime
    time_ago = time() - mtime
    time_str = if time_ago < 3600
        "$(round(Int, time_ago / 60)) minutes ago"
    elseif time_ago < 86400
        "$(round(time_ago / 3600, digits=1)) hours ago"
    else
        "$(round(Int, time_ago / 86400)) days ago"
    end
    
    println("🧠 Found incomplete run ($percent_complete% complete, $(format_number(last_obj))/$(format_number(nobj)) objects)")
    println("   Last updated: $time_str")
    print("   Resume from object $(last_obj + 1)? [Y/n]: ")
    
    response = readline()
    if response == "" || lowercase(response) == "y"
        return ("resume", true)
    else
        return ("overwrite", false)
    end
end

"""
    fit_streaming(param::Dict, templgrid::Array{Float64,3}, template_error_grid::Matrix{Float64},
                  zgrid::Vector{Float64}, bands::Vector{String}, templates::Vector{String},
                  nobj::Int, nz::Int, nband::Int, ntempl::Int, 
                  fnu::Matrix{Float64}, efnu::Matrix{Float64}, IDs::Vector{Int},
                  nphot_min::Int, output_file::String, output_pz::Bool, output_templ::Bool,
                  target_memory_gb::Float64, preserve_work_file::Bool, chunked_processing::Bool)

Main streaming fit function that handles both chunked and in-memory processing modes.
When chunked_processing=false, uses single chunk (in-memory mode).
When chunked_processing=true, uses memory-controlled chunking.
"""
function fit_streaming(param::Dict, templgrid::Array{Float64,3}, template_error_grid::Matrix{Float64},
                      zgrid::Vector{Float64}, bands::Vector{String}, templates::Vector{String},
                      nobj::Int, nz::Int, nband::Int, ntempl::Int, 
                      fnu::Matrix{Float64}, efnu::Matrix{Float64}, IDs::Vector{Int},
                      nphot_min::Int, output_file::String, output_pz::Bool, output_templ::Bool,
                      target_memory_gb::Float64, preserve_work_file::Bool, chunked_processing::Bool)
    
    # Generate work file name
    work_file = output_file * ".work.h5"
    
    # Check for resume capability
    can_resume, last_obj = check_resume_file(work_file)
    start_obj = 1
    
    if can_resume
        action, should_resume = prompt_resume(work_file, last_obj, nobj)
        
        if action == "resume" && should_resume
            start_obj = last_obj + 1
            println("✅ Resuming from object $start_obj")
        elseif action == "complete"
            println("✅ Marking run as complete and generating output...")
            finalize_hdf5_work_file(work_file)
            
            # Jump to output conversion
            final_output = output_file
            if endswith(output_file, ".fits")
                println("💾 Converting HDF5 to FITS format: $output_file")
                convert_hdf5_to_fits(work_file, output_file)
                
                # Handle work file based on preserve_work_file setting
                if preserve_work_file
                    println("💾 Work file preserved: $work_file")
                else
                    print("🗑️  Remove work file $work_file? [Y/n]: ")
                    response = readline()
                    if response == "" || lowercase(response) == "y"
                        rm(work_file)
                        println("✅ Work file removed")
                    else
                        println("💾 Work file preserved: $work_file")
                    end
                end
            else
                # Keep HDF5 format
                if endswith(output_file, ".h5") || endswith(output_file, ".hdf5")
                    mv(work_file, final_output, force=true)
                    println("✅ Results saved to: $final_output")
                else
                    println("✅ Results saved to: $work_file")
                end
            end
            return 0  # Exit early for complete runs
        elseif action == "keep"
            println("💾 Keeping work file and exiting...")
            return 0  # Exit without processing
        else  # action == "overwrite"
            println("🔄 Starting fresh run")
            rm(work_file, force=true)
            can_resume = false
        end
    end
    
    # Create work file if not resuming
    if !can_resume
        println("📝 Creating HDF5 work file: $work_file")
        create_hdf5_work_file(work_file, param, nobj, nz, nband, ntempl, bands, zgrid, templates)
    end
    
    # Determine chunk size based on processing mode
    if chunked_processing
        # Use memory-controlled chunking
        chunk_size = get_chunk_size_for_memory_target(nobj, nz, target_memory_gb)
        actual_memory_gb = chunk_size * nz * 4 / (1024^3)  # Float32 chi2 values
        println("⚙️ Chunked processing: $(format_number(chunk_size)) objects per chunk (~$(round(actual_memory_gb, digits=2)) GB per chunk)")
    else
        # In-memory mode: single chunk = entire catalog
        chunk_size = nobj
        estimated_gb = nobj * nz * 4 / (1024^3)  # Rough estimate for display
        println("⚙️ In-memory processing: all $(format_number(nobj)) objects in single batch (~$(round(estimated_gb, digits=2)) GB)")
        if estimated_gb > 8.0
            println("   ⚠️  High memory usage detected. Consider setting chunked_processing = true if you encounter memory issues")
        end
    end
    
    # Initialize progress tracking
    total_remaining = nobj - start_obj + 1
    progress_bar = ProgressBar(total_remaining, "Fitting objects...", 
                              show_rate=true, show_eta=true, prefix_emoji="🧠")
    
    # Display initial progress at 0%
    display_progress(progress_bar)
    
    # Process objects in chunks
    objects_processed = 0
    
    # Set up interrupt handling
    interrupted = Ref(false)
    
    try
        for chunk_start in start_obj:chunk_size:nobj
            chunk_end = min(chunk_start + chunk_size - 1, nobj)
            chunk_nobj = chunk_end - chunk_start + 1
            
            # Process this chunk
            chunk_IDs = IDs[chunk_start:chunk_end]
            chunk_fnu = fnu[chunk_start:chunk_end, :]
            chunk_efnu = efnu[chunk_start:chunk_end, :]
            
            # Initialize chunk result arrays
            chunk_zbest = zeros(chunk_nobj)
            chunk_chi2best = zeros(chunk_nobj)
            chunk_coeffsbest = zeros(chunk_nobj, ntempl)
            chunk_chi2grid = output_pz ? zeros(Float32, chunk_nobj, nz) : nothing
            
            # Fit objects in this chunk using threading
            nthreads_to_use = min(Threads.nthreads(), chunk_nobj)
            objects_per_task = max(1, chunk_nobj ÷ nthreads_to_use)
            ntasks = cld(chunk_nobj, objects_per_task)
            
            tasks = Vector{Task}(undef, ntasks)
            for task_id in 1:ntasks
                task_start = (task_id - 1) * objects_per_task + 1
                task_end = min(task_id * objects_per_task, chunk_nobj)
                
                tasks[task_id] = Threads.@spawn begin
                    batch_results = Vector{Tuple{Int, Float64, Float64, Vector{Float64}, Vector{Float32}}}()
                    
                    for j_local in task_start:task_end
                        # Check for interruption before starting work
                        if interrupted[]
                            push!(batch_results, (j_local, -1.0, -1.0, zeros(ntempl), fill(Float32(-1.0), nz)))
                            increment!(progress_bar)
                            continue
                        end
                        
                        j_global = chunk_start + j_local - 1
                        
                        zbest_j, chi2best_j, coeffsbest_j, chi2_row_j = fit_single_object(
                            j_global, chunk_fnu[j_local,:], chunk_efnu[j_local,:], 
                            templgrid, template_error_grid, zgrid, nphot_min, 
                            nband, ntempl, nz; interrupted_flag=interrupted
                        )
                        
                        push!(batch_results, (j_local, zbest_j, chi2best_j, coeffsbest_j, chi2_row_j))
                        increment!(progress_bar)
                    end
                    
                    return batch_results
                end
            end
            
            # Collect results from tasks
            for task in tasks
                batch_results = fetch(task)
                for (j_local, zbest_j, chi2best_j, coeffsbest_j, chi2_row_j) in batch_results
                    chunk_zbest[j_local] = zbest_j
                    chunk_chi2best[j_local] = chi2best_j
                    chunk_coeffsbest[j_local, :] = coeffsbest_j
                    if output_pz
                        chunk_chi2grid[j_local, :] = chi2_row_j
                    end
                end
            end
            
            # Calculate P(z) statistics for this chunk
            if output_pz
                pz_chunk = exp.(-0.5 * chunk_chi2grid)
                cpz_chunk = cumsum(pz_chunk, dims=2) ./ sum(pz_chunk, dims=2)
                
                chunk_z_l95 = zgrid[map(argmin, eachrow(abs.(cpz_chunk .- 0.025)))]
                chunk_z_l68 = zgrid[map(argmin, eachrow(abs.(cpz_chunk .- 0.160)))]
                chunk_z_med = zgrid[map(argmin, eachrow(abs.(cpz_chunk .- 0.500)))]
                chunk_z_u68 = zgrid[map(argmin, eachrow(abs.(cpz_chunk .- 0.840)))]
                chunk_z_u95 = zgrid[map(argmin, eachrow(abs.(cpz_chunk .- 0.975)))]
            else
                # Use best-fit redshift as placeholder
                chunk_z_l95 = chunk_zbest
                chunk_z_l68 = chunk_zbest
                chunk_z_med = chunk_zbest
                chunk_z_u68 = chunk_zbest
                chunk_z_u95 = chunk_zbest
                chunk_chi2grid = nothing
            end
            
            # Calculate best-fit photometry for this chunk
            chunk_photobest = zeros(chunk_nobj, nband)
            for j_local in 1:chunk_nobj
                if chunk_zbest[j_local] > 0  # Only for valid fits
                    z_idx = argmin(abs.(zgrid .- chunk_zbest[j_local]))
                    for k in 1:nband
                        chunk_photobest[j_local, k] = 0.0
                        for t in 1:ntempl
                            chunk_photobest[j_local, k] += templgrid[t, z_idx, k] * chunk_coeffsbest[j_local, t]
                        end
                    end
                end
            end
            
            # Write chunk results to HDF5
            write_chunk_results(work_file, chunk_start, chunk_end,
                              chunk_IDs, chunk_zbest, chunk_chi2best, chunk_coeffsbest,
                              chunk_z_l95, chunk_z_l68, chunk_z_med, chunk_z_u68, chunk_z_u95,
                              chunk_photobest, bands, chunk_chi2grid)
            
            # Progress is now updated per-object within the threaded tasks
            objects_processed += chunk_nobj
        end
        
        # Finish progress bar immediately after fitting loop completes
        if interrupted[]
            finish!(progress_bar, final_message="interrupted")
        else
            finish!(progress_bar)
        end
        
        # Generate template output if requested (after fitting is complete)
        if output_templ
            println("📊 Generating template output...")
            
            # Use unified template grid builder for spectral mode
            templgrid_spectral, _, wavelength_grid = build_template_grid(
                templates, zgrid, param["fitting"]["igm_model"], param["fitting"];
                output_type=:spectral
            )
            
            # Write template data to HDF5 work file
            write_template_data(work_file, templgrid_spectral, wavelength_grid, zgrid, templates)
            println("✅ Template output generated")
        end
        
        # Mark work file as complete
        finalize_hdf5_work_file(work_file)
        
        # Handle output format and work file preservation
        final_output = output_file
        if endswith(output_file, ".fits")
            println("💾 Exporting to FITS: $output_file")
            convert_hdf5_to_fits(work_file, output_file)
            
            # Handle work file based on preserve_work_file setting
            if preserve_work_file
                println("💾 Work file preserved: $work_file")
            else
                print("🗑️  Remove work file $work_file? [Y/n]: ")
                response = readline()
                if response == "" || lowercase(response) == "y"
                    rm(work_file)
                    println("✅ Work file removed")
                else
                    println("💾 Work file preserved: $work_file")
                end
            end
        else
            # Keep HDF5 format
            if endswith(output_file, ".h5") || endswith(output_file, ".hdf5")
                mv(work_file, final_output, force=true)
                println("✅ Results saved to: $final_output")
            else
                # Output file doesn't have recognized extension, ask user
                println("⚠️  Output file '$output_file' doesn't have .fits/.h5/.hdf5 extension")
                print("   Save as HDF5 format? [Y/n]: ")
                response = readline()
                if response == "" || lowercase(response) == "y"
                    hdf5_output = output_file * ".h5"
                    mv(work_file, hdf5_output, force=true)
                    println("✅ Results saved to: $hdf5_output")
                else
                    println("💾 Work file preserved: $work_file")
                end
            end
        end
        
    catch e
        if isa(e, InterruptException)
            interrupted[] = true
            println("\n⚠️ Fitting interrupted by user. Processing results for completed objects...")
            finalize_hdf5_work_file(work_file)  # Mark as complete for resume
            finish!(progress_bar, final_message="interrupted")
        else
            println("\n⚠️ Error during streaming fit: $e")
            println("💾 Partial results preserved in: $work_file")
            rethrow(e)
        end
    end
    
    return 0
end