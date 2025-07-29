# Lazy.jl

**Fast, scalable photometric redshift fitting in Julia**

Lazy.jl is a multithreaded photometric redshift fitting package designed for modern astronomical surveys. Built from the ground up in Julia, it provides memory-efficient and multithreaded processing of large catalogs.

[![Julia 1.12+](https://img.shields.io/badge/julia-1.12+-blue.svg)](https://julialang.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Table of Contents

- [Key Features](#key-features)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [CLI Reference](#cli-reference)
- [Configuration](#configuration)
- [Advanced Features](#advanced-features)
- [Output Formats](#output-formats)
- [Performance & Scaling](#performance--scaling)
- [Examples](#examples)
- [Development](#development)

## Key Features

**Performance:**
- Native multithreading using Julia's threading model
- Memory-efficient chunked processing for datasets larger than available RAM
- Template grid caching to reuse computations across runs
- Resume capability for interrupted jobs

**Usability:**
- Simple command-line interface 
- Human-readable TOML configuration files 
- Flexible output formats (FITS or HDF5)
- Real-time progress tracking with memory estimates

**Scientific Features:**
- Built-in intergalactic medium modeling (Inoue+2014)
- Configurable template systematic uncertainties
- Full redshift probability distributions with confidence intervals
- Complete template reconstruction with coefficients

**Customization:**
- Easy addition of custom template sets
- Comprehensive filter database organized by telescope/instrument

## Installation

Lazy.jl requires Julia 1.12+ and is installed from source using the included installation script.

### Prerequisites

1. **Install Julia 1.12+** using [juliaup](https://github.com/JuliaLang/juliaup):
   ```bash
   curl -fsSL https://install.julialang.org | sh
   juliaup add 1.12
   juliaup default 1.12
   ```

2. **Clone and install Lazy.jl**:
   ```bash
   git clone https://github.com/hollisakins/Lazy.jl.git
   cd Lazy.jl
   bash install.sh
   ```

3. **Add to PATH** (add to your `~/.bashrc` or `~/.zshrc`):
   ```bash
   export PATH="$PATH:$HOME/.julia/bin"
   ```

4. **Reload shell**: `source ~/.bashrc` or `source ~/.zshrc`

### Verification

Test your installation:
```bash
lazy --version
lazy list-templates
```

The first run will precompile packages (1-2 minutes), subsequent runs are instant.

### Updates

To update your installation:
```bash
cd Lazy.jl
git pull
bash install.sh
```

## Quick Start

### 1. Create a parameter file

Generate an example parameter file:
```bash
lazy params > my_params.toml
```

### 2. Configure your parameter file

Edit `my_params.toml` to point to your catalog:
```toml
[io]
    input_catalog = '/path/to/your/catalog.fits'
    output_file = '/path/to/results.fits'
    output_pz = true

[fitting]
    template_set = 'sfhz'
    z_min = 0.0
    z_max = 6.0
    z_step = 0.01

[translate]
    # Map Lazy.jl filter names (or nicknames!) to catalog columns
    f606w = {flux = 'f_f606w', error = 'e_f606w'}
    f814w = {flux = 'f_f814w', error = 'e_f814w'}
    # ... add your bands
```

### 3. Run the fit

```bash
lazy fit -p my_params.toml
```

## CLI Reference

Lazy.jl provides a comprehensive command-line interface:

### Main Commands

```bash
lazy fit -p <param_file>          # Run photometric redshift fitting
lazy list-templates               # Show available template sets  
lazy list-filters                 # Show available filter curves
lazy params                       # Export example parameter file
lazy cache-clear                  # Clear cached template grids
lazy --version                    # Show version information
lazy --help                       # Show help information
```

### Fitting Options

```bash
lazy fit -p params.toml           # Basic fitting
lazy fit -p params.toml -t 8      # Use 8 threads
lazy fit -p params.toml -t auto   # Use all available threads
```

### Advanced Usage

```bash
# Resume interrupted job (automatic detection)
lazy fit -p params.toml  # Will prompt to resume if work file exists
```

## Configuration

Lazy.jl uses TOML configuration files organized into logical sections:

### Core Sections

#### `[io]` - Input/Output Configuration
```toml
[io]
    input_catalog = '/path/to/catalog.fits'      # Input FITS catalog
    output_file = '/path/to/results.fits'        # Output file (.fits/.h5/.hdf5)
    missing_data_format = 'nan'                  # Format for missing data
    output_pz = true                             # Include P(z) distributions
    output_templates = false                     # Include best-fit templates
```

#### `[fitting]` - Scientific Parameters
```toml
[fitting]
    template_set = 'sfhz'                        # Template set name or file list
    template_error = 'template_error'            # Template uncertainty model
    template_error_scale = 0.2                   # Template error scaling
    igm_model = 'inoue14'                        # IGM attenuation model
    nphot_min = 2                                # Minimum detection requirement
    sys_err = 0.05                               # Systematic error floor
    z_min = 0.0                                  # Minimum redshift
    z_max = 20.0                                 # Maximum redshift  
    z_step = 0.01                                # Redshift step size
    template_cache = true                        # Enable template caching
```

#### `[runtime]` - Processing Configuration
```toml
[runtime]
    chunked_processing = false                   # Enable memory-efficient mode
    target_memory_gb = 0.5                       # Memory target per chunk
    preserve_work_file = false                   # Keep intermediate files
```

#### `[translate]` - Column Mapping
```toml
[translate]
    # Map filter names to catalog columns
    f606w = {flux = 'f_f606w', error = 'e_f606w'}
    f814w = {flux = 'f_f814w', error = 'e_f814w'}
    # ... continue for all bands in your catalog
```

### Template Sets

Built-in template sets:
- `sfhz`: Star formation history templates (13 templates)
- `fsps_full`: Full FSPS template set (12 templates) 
- `sfhz_Larson22`: Extended SFH templates (15 templates)

Custom template sets:
```toml
[fitting]
    template_set = [
        'custom/template1.dat',
        'custom/template2.fits',
        # ... list all template files
    ]
```

### Filter Management

View available filters:
```bash
lazy list-filters
```

Filters are organized by telescope/instrument:
- `hst/acs/f606w`, `hst/wfc3/f160w`
- `jwst/nircam/f150w`, `jwst/miri/f770w`  
- `spitzer/irac/ch1`, `subaru/hsc/g`
- And many more...

## Advanced Features

### Memory-Efficient Processing

For large catalogs (>100k objects) or large redshift grids, enable chunked processing to keep memory usage manageable. Processing happens in chunks while maintaining multithreading within each chunk. 

```toml
[runtime]
    chunked_processing = true
    target_memory_gb = 1.0                       # Adjust based on your system
```

### Resume Capability

Interrupted jobs can be resumed automatically from HDF5 work files. 

```bash
# If job was interrupted, next run will prompt:
lazy fit -p params.toml
# → "Found incomplete run (75% complete), resume? [Y/n]"
```

Work files (`.work.h5`) store progress and enable resumption from any chunk boundary.

### Template Grid Caching

Template grids are automatically cached to speed up repeated runs:

```bash
lazy cache-clear                                # Clear cache if needed
```

Cache benefits:
- **First run**: Builds and caches template grid (~1m for typical setup)
- **Subsequent runs**: Loads from cache (~2s)
- **Automatic invalidation**: Cache updates when templates/parameters change

## Output Formats

### FITS Output (Default)

Multi-extension FITS file with:

**SUMMARY Extension:**
- `ID`: Object identifiers
- `z_best`: Best-fit redshift
- `chi2`: Best-fit χ²
- `z_l68`, `z_u68`: 68% confidence intervals
- `z_l95`, `z_u95`: 95% confidence intervals  
- `z_med`: Median redshift from P(z)
- `f606w`, `f814w`, ...: Best-fit model fluxes
- `coeffs`: Template coefficients

**PZ Extension** (if `output_pz = true`):
- `ID`: Object identifiers (-1 for the first row, which stores the redshift grid)
- `Pz`: Redshift probability distributions (redshift grid in first row)

**TEMPL Extension** (if `output_templates = true`):
- `z`: Redshift grid (-1 for wavelength grid)
- Template names: Fluxes in each template at each redshift (common wavelength grid in first row)

The `TEMPL` extension allows you to reconstruct the best-fit SED from the `coeffs` for a given object. For example, in python: 

```python
from astropy.io import fits
import numpy as np

ID = 1
lazy_file = '...'

lazy_summary = fits.getdata(lazy_file, ext=1)
i = np.where(lazy_summary['ID']==ID)[0][0]
z_best = lazy_summary['z_best'][i]
coeffs = lazy_summary['coeffs'][i]

lazy_templ = fits.getdata(lazy_file, ext=3)
izbest = np.argmin(np.abs(z_best-lazy_templ['z']))
template_names = [n for n in lazy_templ.dtype.names if n!='z']
templates = np.array([lazy_templ[template][izbest] for template in template_names])

wavelength = lazy_templ[template_names[0]][0] * (1+z_best) # Wavelength in first row
fnu = np.dot(templates.T, coeffs) 
```

### HDF5 Output

For streaming mode or when output file ends in `.h5`/`.hdf5`:

```
results.h5
├── metadata/           # Run parameters and provenance
├── results/            # Main fitting results  
├── photometry/         # Best-fit model photometry
└── pz/                 # P(z) distributions (compressed)
```

## Performance & Scaling

### Computational Complexity

- **Template grid construction**: O(N_templates × N_redshifts × N_bands)
- **Object fitting**: O(N_objects × N_redshifts × N_templates × N_bands)
- **Memory usage**: O(N_objects × N_redshifts) for P(z) storage

### Threading Performance

Scaling with thread count (1M objects, 10 bands, 20 templates):

| Threads | Time    | Speedup | Efficiency |
|---------|---------|---------|------------|
| 1       | 45 min  | 1.0x    | 100%       |
| 4       | 12 min  | 3.8x    | 95%        |
| 8       | 6.5 min | 6.9x    | 86%        |
| 16      | 4.2 min | 10.7x   | 67%        |

### Memory Usage Examples

| Dataset Size | Objects | Redshifts | Memory (In-memory) | Memory (Chunked) |
|--------------|---------|-----------|-------------------|------------------|
| Small        | 10k     | 200       | 0.02 GB          | N/A              |
| Medium       | 100k    | 600       | 0.45 GB          | N/A              |
| Large        | 500k    | 2000      | 3.7 GB           | 1.0 GB           |
| Very Large   | 1M      | 2000      | 7.5 GB           | 1.0 GB           |
| Extreme      | 10M     | 2000      | 75 GB            | 1.0 GB           |

### Optimization Tips

1. **Choose appropriate z_step**: Fine grids (0.001) vs coarse grids (0.01) have 100x memory difference
2. **Use chunked processing**: For datasets >1GB estimated memory usage
3. **Template caching**: Significant speedup for repeated runs with same templates
4. **Thread tuning**: Optimal thread count is usually 0.5-1x CPU cores
5. **Memory targets**: Match `target_memory_gb` to available system RAM

## Examples

### Basic Photometric Redshift Fitting

```toml
# basic_fitting.toml
[io]
    input_catalog = 'cosmos_catalog.fits'
    output_file = 'cosmos_photoz.fits'
    output_pz = true
    output_templates = true

[fitting]
    template_set = 'sfhz'
    template_error = 'template_error'
    template_error_scale = 0.2
    igm_model = 'inoue14'
    nphot_min = 2 
    z_min = 0.0
    z_max = 6.0
    z_step = 0.03                               # Coarse grid for speed
    sys_err = 0.05

[translate]
    f606w = {flux = 'f_auto_f606w', error = 'e_auto_f606w'}
    f814w = {flux = 'f_auto_f814w', error = 'e_auto_f814w'}
    f125w = {flux = 'f_auto_f125w', error = 'e_auto_f125w'}
    f160w = {flux = 'f_auto_f160w', error = 'e_auto_f160w'}
```

```bash
lazy fit -p basic_fitting.toml
```

### Getting Help

1. **Check parameter file**: Run `lazy params` for working example
2. **Verify installation**: Run `lazy --version` and `lazy list-templates`
3. **Check logs**: Look for warning messages in output
4. **GitHub Issues**: Report bugs at https://github.com/hollisakins/Lazy.jl/issues

## Development

### Contributing

Contributions are welcome! Please see our contribution guidelines:

1. **Fork** the repository
2. **Create** a feature branch: `git checkout -b feature-name`
3. **Test** your changes: `julia --project=. -e "using Pkg; Pkg.test()"`
4. **Submit** a pull request

### Architecture Overview

```
src/
├── Lazy.jl              # Main module definition and version checking
├── cli.jl               # Command-line interface and user interaction
├── fitting.jl           # Core photometric redshift fitting algorithms
├── io.jl                # I/O operations, caching, and data management
├── utils.jl             # Utilities (progress bars, memory estimation, formatting)
├── template_grid.jl     # Template grid construction and IGM modeling
├── writedata.py         # Python interface for FITS I/O
├── templates/           # SED template library
├── filter_files/        # Filter transmission curves
└── igm_data/           # IGM attenuation models
```

### Adding Templates

1. **Copy template files** to `src/templates/your_set/`
2. **Update template directory**:
   ```toml
   # src/templates/template_directory.toml
   [your_set]
   files = ['your_set/template1.dat', 'your_set/template2.dat']
   description = 'Your custom template set'
   ```
3. **Test**: `lazy list-templates` should show your set

### Adding Filters

1. **Copy filter file** to `src/filter_files/telescope/instrument/`
2. **Update filter directory**:
   ```toml  
   # src/filter_files/filter_directory.toml
   ['telescope/instrument/filter']
   description = 'Description of filter'
   nicknames = ['shortname', 'alt_name']
   ```
3. **Test**: `lazy list-filters` should show your filter

### Extending Functionality

The modular architecture makes it easy to extend Lazy.jl:

- **New IGM models**: Add to `src/igm_data/`
- **Template error functions**: Add to `src/templates/template_error/`
- **Output formats**: Extend I/O functions in `io.jl`
- **Fitting algorithms**: Modify core algorithms in `fitting.jl`

### Testing

Run the test suite:
```bash
julia --project=. -e "using Pkg; Pkg.test()"
```

Performance benchmarks:
```bash
julia --project=. test_phase1_demo.jl
julia --project=. runtime_config_demo.jl
```

---

## License

Lazy.jl is released under the MIT License. See [LICENSE](LICENSE) for details.

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

## Acknowledgments

- Based on the [eazy-py](https://github.com/gbrammer/eazy-py) package by Gabriel Brammer
- IGM attenuation from [Inoue et al. 2014](https://ui.adsabs.harvard.edu/abs/2014MNRAS.442.1805I)
- Built with the [Julia](https://julialang.org) programming language

---

*Lazy.jl: Making photometric redshift fitting fast, scalable, and user-friendly.*