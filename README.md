# Lazy.jl



[![Build Status](https://github.com/hollisakins/Lazy.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/hollisakins/Lazy.jl/actions/workflows/CI.yml?query=branch%3Amain)




## Installation 

Still trying to find a clean way to install as a command line executable. For now, the best solution I've come up with: 

Make sure you have `julia` installed. 

### Install Lazy.jl
```
git clone https://github.com/hollisakins/Lazy.jl.git
cd Lazy.jl
julia --project=@lazy -e 'using Pkg; Pkg.develop(path=".")'
```
### Install the lazy shell script
```
install="~/.julia/bin/lazy"
curl -fsSL -o $install https://raw.githubusercontent.com/.../.../lazy
chmod +x $install
```
