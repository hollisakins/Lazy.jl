module Lazy

export LazyError, main, templatepath, filterpath

const version = string(pkgversion(@__MODULE__))

using LinearAlgebra: BLAS
BLAS.set_num_threads(1)

using TOML, HTTP
function check_version()
    # Check if we need to do a version check based on last check time
    lazy_dir = joinpath(homedir(), ".lazy")
    version_file = joinpath(lazy_dir, "version")
    
    # Only check version if cache file doesn't exist or is older than 1 day
    if isfile(version_file)
        cache_time = mtime(version_file)
        if time() - cache_time < 24 * 3600  # 1 day in seconds
            return 0  # Skip check, too recent
        end
    end
    
    # Perform the actual version check
    url = "https://raw.githubusercontent.com/hollisakins/Lazy.jl/main/Project.toml"
    try
        r = HTTP.get(url)
        if r.status == 200
            latest_version = TOML.parse(String(r.body))["version"]
            if version != latest_version
                @warn "Lazy.jl is out of date. Please update to version $latest_version (e.g. `git pull` and `bash install.sh`)."
                # Update version file even if out of date
                mkpath(lazy_dir)
                touch(version_file)
                return 1
            else
                # Update version file to record successful check
                mkpath(lazy_dir)
                touch(version_file)
            end
        end
    catch e
        # Silently ignore network errors, don't update cache
        # Could optionally log debug info: @debug "Version check failed: $e"
    end
    return 0
end


using CFITSIO

include("template_grid.jl")
include("utils.jl")
include("io.jl")
include("fitting.jl")
include("cli.jl")

# Precompile the entry points
@assert precompile(main, (Vector{String},))

end