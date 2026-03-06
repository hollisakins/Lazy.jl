# Lazy.jl

**Fast, scalable photometric redshift fitting in Julia**

Lazy.jl is a multithreaded photometric redshift fitting package designed for modern astronomical surveys. Built from the ground up in Julia, it provides memory-efficient and multithreaded processing of large catalogs using NNLS template fitting of SED models with IGM attenuation.

## Why Lazy.jl?

Lazy.jl is a Julia reimplementation of [eazy-py](https://github.com/gbrammer/eazy-py) focused on performance and ease of use.

- **Fast**: Native multithreading scales near-linearly across cores, with shared memory instead of the per-process overhead of Python multiprocessing.
- **Scalable**: Chunked processing fits arbitrarily large catalogs within a fixed memory budget, with automatic resume for interrupted jobs.
- **Simple**: A single CLI (`lazy fit -p params.toml`) with human-readable TOML configuration, easy to integrate into existing pipelines.

## Quick Start

```bash
# Install
git clone https://github.com/hollisakins/Lazy.jl.git
cd Lazy.jl && bash install.sh

# Generate a parameter file and run
lazy params > my_params.toml
# Edit my_params.toml to point to your catalog...
lazy fit -p my_params.toml
```

See [Installation](@ref) for detailed setup instructions and [Getting Started](@ref) for a step-by-step tutorial.
