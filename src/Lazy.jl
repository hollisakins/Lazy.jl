module Lazy

const version = pkgversion(@__MODULE__)

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