# Lazy.jl

**Fast, scalable photometric redshift fitting in Julia**

Lazy.jl is a multithreaded photometric redshift fitting package designed for modern astronomical surveys. Built from the ground up in Julia, it provides memory-efficient processing of large catalogs with built-in fault tolerance and resume capability.

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
- [Technical Details](#technical-details)
- [Troubleshooting](#troubleshooting)
- [Development](#development)

## Key Features

### 🚀 **High Performance**
- **Native multithreading**: Shared memory parallelization using Julia's threading model
- **Memory-efficient streaming**: Process datasets larger than available RAM
- **Template grid caching**: Reuse computed grids across runs
- **Optimized algorithms**: Fast non-negative least squares fitting with numerical stability

### 📊 **Scalable Architecture**
- **Chunked processing**: Automatic memory management for large catalogs (1M+ objects)
- **Resume capability**: Interrupted jobs can be resumed from checkpoints
- **Memory targets**: User-configurable memory usage (0.1GB to 10GB+ per chunk)
- **Progress tracking**: Real-time progress with ETA estimation

### 🔧 **User-Friendly**
- **Simple CLI**: Intuitive command-line interface with comprehensive help
- **TOML configuration**: Human-readable parameter files with clear documentation
- **Flexible I/O**: FITS and HDF5 output formats with automatic conversion
- **Rich feedback**: Detailed logging, memory estimates, and performance metrics

### 🧬 **Scientific Accuracy**
- **IGM attenuation**: Built-in intergalactic medium modeling (Inoue+2014)
- **Template errors**: Configurable systematic uncertainties
- **P(z) distributions**: Full redshift probability distributions with confidence intervals
- **Best-fit SEDs**: Complete template reconstruction with coefficients

### 🏗️ **Extensible Design**
- **Template support**: Easy addition of custom SED template sets
- **Filter management**: Comprehensive filter database with telescope/instrument organization
- **Modular architecture**: Clean separation of fitting, I/O, and processing components

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

### 2. Configure your data

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
    # Map filter names to catalog columns
    f606w = {flux = 'f_f606w', error = 'e_f606w'}
    f814w = {flux = 'f_f814w', error = 'e_f814w'}
    # ... add your bands
```

### 3. Run the fit

```bash
lazy fit -p my_params.toml
```

For large datasets, enable memory-efficient processing:
```toml
[runtime]
    chunked_processing = true
    target_memory_gb = 1.0
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
# Check memory requirements before running
lazy fit -p params.toml --dry-run

# Verbose output for debugging
lazy fit -p params.toml --verbose

# Resume interrupted job (automatic if work file exists)
lazy fit -p params.toml  # Will prompt to resume if applicable
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

For large catalogs (>100k objects), enable chunked processing:

```toml
[runtime]
    chunked_processing = true
    target_memory_gb = 1.0                       # Adjust based on your system
```

**Memory targets guide:**
- `0.25 GB`: Conservative (many small chunks)
- `0.5 GB`: Balanced (recommended default)
- `1.0 GB`: Moderate (good performance)
- `2.0+ GB`: Aggressive (requires sufficient RAM)

### Resume Capability

Interrupted jobs can be resumed automatically:

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
- **First run**: Builds and caches template grid (~30s for typical setup)
- **Subsequent runs**: Loads from cache (~2s)
- **Automatic invalidation**: Cache updates when templates/parameters change

### Performance Monitoring

Lazy.jl provides detailed performance feedback:

```
📊 Objects: 1.0M
📊 Bands: 10  
📊 Filters: f606w, f814w, f090w, f115w, f150w, f200w, f277w, f356w, f410m, f444w

====================================
Memory Usage Estimate:
  Template grid:      0.0 GB
  Chi2 grid:          7.5 GB
  Template errors:    0.0 GB
  Coefficients:       0.15 GB
  Photometric data:   0.15 GB
  Working arrays:     0.0 GB
  ----------------------------------------
  Estimated peak:     7.8 GB
  ⚠️  CAUTION: Moderate memory usage.
     Monitor system memory during fitting.
====================================

⚙️ Runtime: Chunked processing enabled
   Target memory: 1.0 GB per chunk
⚙️ Using chunks of 134.2k objects (~1.0 GB per chunk)

🧠 Fitting objects... ████████████████████ 100% (1.0M/1.0M; 8.4k obj/s, ETA: 0s)
```

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
- `ID`: Object identifiers (-1 for redshift grid)
- `Pz`: Redshift probability distributions

**TEMPL Extension** (if `output_templates = true`):
- `z`: Redshift grid (-1 for wavelength grid)
- Template names: Best-fit SED templates

### HDF5 Output

For streaming mode or when output file ends in `.h5`/`.hdf5`:

```
results.h5
├── metadata/           # Run parameters and provenance
├── results/            # Main fitting results  
├── photometry/         # Best-fit model photometry
└── pz/                 # P(z) distributions (compressed)
```

HDF5 benefits:
- **Compression**: P(z) data compressed ~3-5x
- **Metadata**: Complete parameter provenance
- **Inspection**: Can be examined with standard HDF5 tools
- **Streaming**: Enables processing of arbitrarily large datasets

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

1. **Choose appropriate z_step**: Fine grids (0.01) vs coarse grids (0.1) have 100x memory difference
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

[fitting]
    template_set = 'sfhz'
    z_min = 0.0
    z_max = 6.0
    z_step = 0.05                               # Coarse grid for speed

[translate]
    f606w = {flux = 'f_auto_f606w', error = 'e_auto_f606w'}
    f814w = {flux = 'f_auto_f814w', error = 'e_auto_f814w'}
    f125w = {flux = 'f_auto_f125w', error = 'e_auto_f125w'}
    f160w = {flux = 'f_auto_f160w', error = 'e_auto_f160w'}
```

```bash
lazy fit -p basic_fitting.toml
```

### High-Precision Fitting for Publication

```toml
# precision_fitting.toml  
[io]
    input_catalog = 'final_catalog.fits'
    output_file = 'publication_photoz.fits'
    output_pz = true
    output_templates = true

[fitting]
    template_set = 'sfhz_Larson22'              # Extended template set
    template_error_scale = 0.15                 # Conservative template errors
    z_min = 0.0
    z_max = 15.0  
    z_step = 0.01                               # Fine redshift grid
    nphot_min = 3                               # Strict detection requirement
    sys_err = 0.03                              # Low systematic error

[runtime]
    chunked_processing = true                   # Required for fine grid
    target_memory_gb = 2.0                      # Use available memory
    preserve_work_file = true                   # Keep for analysis

[translate]
    # Full HST+JWST coverage
    f435w = {flux = 'f_auto_f435w', error = 'e_auto_f435w'}
    f606w = {flux = 'f_auto_f606w', error = 'e_auto_f606w'}
    f775w = {flux = 'f_auto_f775w', error = 'e_auto_f775w'}
    f814w = {flux = 'f_auto_f814w', error = 'e_auto_f814w'}
    f850lp = {flux = 'f_auto_f850lp', error = 'e_auto_f850lp'}
    f090w = {flux = 'f_auto_f090w', error = 'e_auto_f090w'}
    f115w = {flux = 'f_auto_f115w', error = 'e_auto_f115w'}
    f150w = {flux = 'f_auto_f150w', error = 'e_auto_f150w'}
    f200w = {flux = 'f_auto_f200w', error = 'e_auto_f200w'}
    f277w = {flux = 'f_auto_f277w', error = 'e_auto_f277w'}
    f356w = {flux = 'f_auto_f356w', error = 'e_auto_f356w'}
    f444w = {flux = 'f_auto_f444w', error = 'e_auto_f444w'}
```

### Large Survey Processing

```toml
# survey_processing.toml
[io]
    input_catalog = 'survey_10million.fits'
    output_file = 'survey_photoz.h5'            # Use HDF5 for large datasets
    output_pz = false                           # Skip P(z) to save space/time

[fitting]
    template_set = 'sfhz'
    z_min = 0.0
    z_max = 8.0
    z_step = 0.02                               # Balanced resolution
    nphot_min = 2                               # Inclusive for large surveys

[runtime]
    chunked_processing = true
    target_memory_gb = 4.0                      # Use available memory
    preserve_work_file = false                  # Auto-cleanup
```

## Technical Details

### SED Fitting Algorithm

Lazy.jl implements χ² minimization with non-negative least squares:

1. **Template Grid Construction**: 
   - Load SED templates and apply IGM attenuation (Inoue+2014)
   - Integrate through filter transmission curves
   - Cache results for reuse

2. **Object Fitting**:
   - For each object and redshift, solve: `min ||A·c - b||²` where `c ≥ 0`
   - `A`: template fluxes, `b`: observed fluxes, `c`: coefficients
   - Use numerically stable NNLS algorithm

3. **Statistical Analysis**:
   - Convert χ² to P(z): `P(z) ∝ exp(-χ²/2)`
   - Compute confidence intervals from cumulative P(z)
   - Calculate best-fit model photometry

### Numerical Stability

- **Template normalization**: Templates normalized to rest-frame 1μm
- **Error handling**: Robust handling of non-finite values and edge cases
- **Precision**: Float32 for P(z) storage, Float64 for computations
- **Validation**: Comprehensive input validation and error reporting

### Memory Architecture

- **Shared template grid**: Read-only access across all threads
- **Per-thread working arrays**: Avoid memory allocation in hot loops
- **Chunked P(z)**: Process and write P(z) in manageable chunks
- **Streaming I/O**: HDF5 backend for memory-efficient operations

## Troubleshooting

### Common Issues

**Memory Errors**
```
ERROR: OutOfMemoryError()
```
**Solution**: Enable chunked processing or reduce `target_memory_gb`

**Template Not Found**
```
LazyError: Template set 'custom' not found in template directory
```
**Solution**: Check template names with `lazy list-templates` or verify file paths

**Column Missing**
```
LazyError: Column 'f_f606w' not found in catalog
```
**Solution**: Verify column names in your FITS catalog and update `[translate]` section

**No Valid Photometry**
```
Warning: No valid bands found for object 12345
```
**Solution**: Check `missing_data_format` and `nphot_min` settings

### Performance Issues

**Slow Startup**
- First run requires package precompilation (1-2 minutes)
- Subsequent runs are fast
- Use `julia --project=. -e "using Lazy"` to precompile manually

**High Memory Usage**
- Check memory estimate before running
- Use `target_memory_gb` to control chunk size
- Consider coarser redshift grid for initial exploration

**Poor Threading Performance**
- Optimal thread count is usually 0.5-1x CPU cores
- BLAS threads are automatically set to 1 (avoid oversubscription)
- Very small datasets may not benefit from many threads

### Debugging

**Enable verbose output**:
```bash
export JULIA_DEBUG=Lazy
lazy fit -p params.toml
```

**Check work files**:
```bash
# Examine partial results
h5dump results.work.h5
```

**Memory monitoring**:
```bash
# Linux/macOS
htop
# Monitor memory usage during fitting
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
├── Lazy.jl              # Main module and CLI entry point
├── main.jl              # Core fitting algorithms and CLI parsing
├── hdf5_streaming.jl    # Streaming processing and HDF5 I/O
├── cache_utils.jl       # Template grid caching system
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
- **Output formats**: Extend I/O functions in `hdf5_streaming.jl`
- **Fitting algorithms**: Modify core algorithms in `main.jl`

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