module Lazy

include("main.jl")

# Precompile the entry points
@assert precompile(main, (Cint,))

end