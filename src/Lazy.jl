module Lazy

export LazyError, main, templatepath, filterpath

const version = string(pkgversion(@__MODULE__))

using LinearAlgebra: BLAS
BLAS.set_num_threads(1)

using TOML, HTTP
function check_version()
    # some logic to ensure the code is up-to-date
    url = "https://raw.githubusercontent.com/hollisakins/Lazy.jl/main/Project.toml"
    r = HTTP.get(url)
    if r.status == 200
        latest_version = TOML.parse(String(r.body))["version"]
        if version != latest_version
            println("Warning: Lazy.jl is out of date. Please update to version $latest_version (e.g. `git pull` and `bash install.sh`).")
            return 1
        end
    end
    return 0
end


using PyCall
const writedata = PyNULL()
function __init__()
    pyimport_conda("astropy.io", "astropy")
    pushfirst!(pyimport("sys")."path", @__DIR__())
    copy!(writedata, pyimport("writedata"))
end

function write_data(filename, data, extname)
    return writedata.write_data(filename, data, extname)
end

include("cache_utils.jl")
include("main.jl")

# CLI entry point with proper exception handling
function julia_main()::Cint
    try
        return main(copy(ARGS))
    catch e
        if isa(e, LazyError)
            printstyled(stderr, "ERROR: "; color = :red, bold = true)
            println(stderr, e.msg)
            return 1
        else
            printstyled(stderr, "UNEXPECTED ERROR: "; color = :red, bold = true)
            println(stderr, sprint(showerror, e))
            Base.show_backtrace(stderr, catch_backtrace())
            return 1
        end
    end
end

# Precompile the entry points
@assert precompile(main, (Vector{String},))
@assert precompile(julia_main, ())

end