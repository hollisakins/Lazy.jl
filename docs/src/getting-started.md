# Getting Started

This guide walks you through your first photometric redshift fit with Lazy.jl.

## 1. Generate a Parameter File

Lazy.jl uses TOML configuration files. Generate a starting template:

```bash
lazy params > my_params.toml
```

This creates a fully commented example file with all available options.

## 2. Configure Your Parameter File

Open `my_params.toml` in a text editor. The file has four sections:

### `[io]` -- Input and Output

Point to your FITS catalog and choose an output location:

```toml
[io]
    input_catalog = '/path/to/your/catalog.fits'
    output_file = '/path/to/results.fits'
    missing_data_format = 'nan'
    output_pz = true
    output_templates = true
```

### `[fitting]` -- Scientific Parameters

Choose a template set and configure the redshift grid:

```toml
[fitting]
    template_set = 'sfhz'
    template_error = 'template_error'
    template_error_scale = 0.2
    igm_model = 'inoue14'
    nphot_min = 2
    sys_err = 0.05
    z_min = 0.0
    z_max = 6.0
    z_step = 0.01
```

Use `lazy list-templates` to see all available template sets.

### `[runtime]` -- Processing Options

For most datasets, the defaults are fine. Enable chunked processing for very large catalogs:

```toml
[runtime]
    chunked_processing = false
    target_memory_gb = 0.5
```

### `[translate]` -- Column Mapping

This is the most important section to get right. Map Lazy.jl filter names to the flux and error columns in your catalog:

```toml
[translate]
    f606w = {flux = 'f_f606w', error = 'e_f606w'}
    f814w = {flux = 'f_f814w', error = 'e_f814w'}
    f090w = {flux = 'f_f090w', error = 'e_f090w'}
    f115w = {flux = 'f_f115w', error = 'e_f115w'}
    f150w = {flux = 'f_f150w', error = 'e_f150w'}
    f200w = {flux = 'f_f200w', error = 'e_f200w'}
    f277w = {flux = 'f_f277w', error = 'e_f277w'}
    f356w = {flux = 'f_f356w', error = 'e_f356w'}
    f410m = {flux = 'f_f410m', error = 'e_f410m'}
    f444w = {flux = 'f_f444w', error = 'e_f444w'}
```

The keys (e.g., `f606w`) are filter nicknames -- use `lazy list-filters` to see all available filters and their nicknames. The values must match column names in your input FITS catalog.

## 3. Run the Fit

```bash
lazy fit -p my_params.toml
```

To use multiple threads for faster processing:

```bash
lazy fit -p my_params.toml -t 8      # Use 8 threads
lazy fit -p my_params.toml -t auto   # Use all available threads
```

Lazy.jl will display a progress bar with estimated time remaining and memory usage.

## 4. Inspect the Output

The output FITS file contains multiple extensions:

- **SUMMARY**: Best-fit redshifts, chi-squared values, confidence intervals, and model photometry for every object
- **PZ** (if `output_pz = true`): Full redshift probability distributions
- **TEMPL** (if `output_templates = true`): Template SEDs with IGM attenuation, enabling SED reconstruction

You can read the results in Python:

```python
from astropy.io import fits

results = fits.getdata('results.fits', ext=1)
print(results['z_best'])    # Best-fit redshifts
print(results['z_med'])     # Median redshifts from P(z)
print(results['z_l68'])     # Lower 68% confidence bound
print(results['z_u68'])     # Upper 68% confidence bound
```

See [Output Formats](@ref) for full details on all output columns and how to reconstruct best-fit SEDs.

## Next Steps

- [Configuration](@ref): Complete reference for all parameter file options
- [Templates](@ref): Learn about the available template sets
- [Filters](@ref): Browse the full filter database
- [Advanced Usage](@ref): Chunked processing, caching, resume, and more
