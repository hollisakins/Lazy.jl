# Lazy.jl

**Fast, scalable photometric redshift fitting in Julia**

[![Julia 1.12+](https://img.shields.io/badge/julia-1.12+-blue.svg)](https://julialang.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://hollisakins.com/Lazy.jl/)

Lazy.jl is a multithreaded photometric redshift fitting package designed for modern astronomical surveys. Built from the ground up in Julia, it provides memory-efficient and multithreaded processing of large catalogs using SED template fitting with IGM attenuation models. Native multithreading scales near-linearly across cores — Lazy can fit up to ~3k obj/s on modern Apple Silicon (16 cores), or up to ~10k obj/s on enterprise chips (64 cores) — churning through catalogs of >1million objects in minutes. 

## Quick Start

```bash
# Install Julia 1.12+ (https://github.com/JuliaLang/juliaup)
# Then:
git clone https://github.com/hollisakins/Lazy.jl.git
cd Lazy.jl && bash install.sh

# Generate a parameter file, edit it, and run
lazy params > my_params.toml
lazy fit -p my_params.toml -t auto
```

## Documentation

Full documentation is available at **[hollisakins.com/Lazy.jl/](https://hollisakins.com/Lazy.jl/)**, including:

- [Installation guide](https://hollisakins.com/Lazy.jl/installation/)
- [Getting started tutorial](https://hollisakins.com/Lazy.jl/getting-started/)
- [Configuration reference](https://hollisakins.com/Lazy.jl/configuration/)
- [CLI reference](https://hollisakins.com/Lazy.jl/cli/)
- [Template and filter catalogs](https://hollisakins.com/Lazy.jl/templates/)
- [Output format details](https://hollisakins.com/Lazy.jl/output/)

## Citation

If you use Lazy.jl in your research, please cite:

```bibtex
@software{lazy_jl,
  title = {Lazy.jl: Fast photometric redshift fitting in Julia},
  author = {Akins, Hollis B.},
  year = {2024},
  url = {https://github.com/hollisakins/Lazy.jl}
}
```

## License

Lazy.jl is released under the MIT License. See [LICENSE](LICENSE) for details.

## Acknowledgments

- Based on the [eazy-py](https://github.com/gbrammer/eazy-py) package by Gabriel Brammer
- IGM attenuation from [Inoue et al. 2014](https://ui.adsabs.harvard.edu/abs/2014MNRAS.442.1805I)
- CGM damping wing model from [Asada et al. 2024](https://ui.adsabs.harvard.edu/abs/2024ApJ...976..140A)
- Built with the [Julia](https://julialang.org) programming language
