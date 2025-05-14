module Lazy

const version = string(pkgversion(@__MODULE__))

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

function write_data(filename, data; extname="DATA")
    return writedata.write_data(filename, data, extname=extname)
end

include("main.jl")

# Precompile the entry points
@assert precompile(main, (Cint,))

end