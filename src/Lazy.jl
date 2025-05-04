module Lazy

const version = pkgversion(@__MODULE__)

include("main.jl")

# Precompile the entry points
@assert precompile(main, (Cint,))

end