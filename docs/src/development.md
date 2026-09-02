# Development

## Architecture

```
src/
├── Lazy.jl              # Main module definition and version checking
├── cli.jl               # Command-line interface and user interaction
├── fitting.jl           # Core photometric redshift fitting algorithms
├── io.jl                # I/O operations, caching, and data management
├── utils.jl             # Utilities (progress bars, memory estimation, formatting)
├── template_grid.jl     # Template grid construction and IGM/CGM modeling
├── example_params.toml  # Example parameter file (used by `lazy params`)
├── templates/           # SED template library
│   └── template_directory.toml
├── filter_files/        # Filter transmission curves
│   └── filter_directory.toml
└── igm_data/            # IGM attenuation model data (Inoue+2014)
```

### Data Flow

1. **Input**: FITS catalog + TOML parameter file
2. **Template grid construction** (`template_grid.jl`): Interpolate templates onto the redshift grid, apply IGM/CGM attenuation, integrate through filter bandpasses
3. **Fitting** (`fitting.jl`): For each object, solve non-negative least squares at every redshift to find the best-fit template combination
4. **Output** (`io.jl`): Write results to FITS (via CFITSIO.jl) or HDF5. Lazy.jl is pure Julia and has no Python dependency.

### Key Algorithms

- **Non-negative least squares** (`NonNegLeastSquares.jl`): Ensures template coefficients are non-negative (physically motivated -- you can't have negative flux contributions)
- **IGM attenuation**: Inoue et al. (2014) model with Lyman-series and continuum absorption
- **CGM damping wing**: Asada et al. (2024) Lyman-alpha damping wing model
- **P(z) computation**: exp(-0.5 * chi2) normalized over the redshift grid

## Contributing

Contributions are welcome, especially new filter transmission curves and template sets! If you've added filters for an instrument or survey not yet included, or have a useful template set, please submit a PR to make them available to the community.

1. **Fork** the repository
2. **Create** a feature branch: `git checkout -b feature-name`
3. **Test** your changes: `julia --project=. -e "using Pkg; Pkg.test()"`
4. **Submit** a pull request

## Running Tests

```bash
julia --project=. -e "using Pkg; Pkg.test()"
```

The test suite covers CLI parsing, template/filter directory integrity, CGM model physics, spectroscopic redshift fixing, P(z) numerics, and the HDF5 to FITS export.

## Extending Lazy.jl

### Adding Templates

See [Templates](@ref) for the full guide. In brief: add files to `src/templates/` and register them in `template_directory.toml`.

### Adding Filters

See [Filters](@ref) for the full guide. In brief: add transmission curves to `src/filter_files/` and register them in `filter_directory.toml`.

### Other Extension Points

- **IGM models**: Add attenuation data to `src/igm_data/` and implement the loader
- **Template error functions**: Add to `src/templates/template_error/`
- **Output formats**: Extend I/O functions in `src/io.jl`
- **Fitting algorithms**: Modify the core fitting loop in `src/fitting.jl`
