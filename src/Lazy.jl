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
function __init__()
    pyimport_conda("astropy.table", "astropy")

    py"""
    from astropy.table import Table
    def write_data(filename, columns, data):
        pydata = {}
        for i, col in enumerate(columns):
            pydata[col] = data[i]
        t = Table(pydata)
        t.write(filename, format='fits', overwrite=True)
        return
    """
end

function write_data(filename, columns, data)
    py"write_data"(filename, columns, data)
end

include("main.jl")

# Precompile the entry points
@assert precompile(main, (Cint,))

end