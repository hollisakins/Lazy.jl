# Configuration

Lazy.jl uses TOML parameter files with four sections. Generate an annotated example with `lazy params > example.toml`.

## `[io]` -- Input/Output

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `input_catalog` | String | *required* | Path to input FITS catalog |
| `output_file` | String | *required* | Output file path (`.fits`, `.h5`, or `.hdf5`) |
| `output_format` | String | inferred | Output format: `'fits'`, `'hdf5'`, or `'both'`. If omitted, inferred from `output_file` extension |
| `missing_data_format` | String | *required* | How missing data is encoded in the catalog: `'nan'`, `'zero'`, or a numeric value (e.g., `'-99'`) |
| `output_pz` | Bool | `false` | Include P(z) distributions in output |
| `output_templates` | Bool | `false` | Include IGM-attenuated template SEDs in output (needed for SED reconstruction) |

### Rest-Frame Magnitudes

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `output_restframe_mags` | Bool | `false` | Compute rest-frame absolute magnitudes (M\_UV, M\_U, M\_V, M\_J) |
| `flux_units` | String | -- | **Required if `output_restframe_mags = true`**. Flux units of the input catalog: `'uJy'`, `'nJy'`, `'Jy'`, or `'cgs'` |
| `H0` | Float | `70.0` | Hubble constant in km/s/Mpc (flat LCDM cosmology) |
| `Om` | Float | `0.3` | Matter density parameter (flat LCDM cosmology) |

### Forced Low-Redshift Fitting

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `output_forced_lowz` | Bool | `false` | Perform a second fit restricted to low redshifts |
| `forced_lowz_zmax` | Float | -- | **Required if `output_forced_lowz = true`**. Maximum redshift for the forced low-z fit |

When enabled, the output includes `z_best_lowz`, `chi2_lowz`, `delta_chi2` (= chi2\_lowz - chi2\_best), low-z P(z) quantiles, and low-z model photometry. This is useful for validating high-redshift candidates.

## `[fitting]` -- Scientific Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `template_set` | String or Array | *required* | Template set name (e.g., `'sfhz'`) or list of custom template file paths |
| `template_error` | String | *required* | Template error model. Currently only `'template_error'` is supported |
| `template_error_scale` | Float | *required* | Scale factor for the template error function |
| `igm_model` | String | *required* | IGM attenuation model. Currently only `'inoue14'` is supported |
| `nphot_min` | Int | *required* | Minimum number of valid photometric points (SNR > 2) required to run the fit |
| `sys_err` | Float | *required* | Systematic error added in quadrature to photometric errors |
| `z_min` | Float | *required* | Minimum redshift for the fitting grid |
| `z_max` | Float | *required* | Maximum redshift for the fitting grid |
| `z_step` | Float | *required* | Redshift step size |
| `template_cache` | Bool | `true` | Cache template grids for reuse in subsequent runs |

### Spectroscopic Redshifts

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `use_zspec` | Bool | `false` | Fix the fit at the spectroscopic redshift when available |

When enabled, objects with a valid spectroscopic redshift (not NaN, not negative) are fit only at the nearest grid point to their z\_spec. Objects without z\_spec are fit normally across the full grid. Requires a `zspec` entry in the `[translate]` section.

### CGM Damping Wing

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `add_cgm` | Bool | `true` | Enable CGM Lyman-alpha damping wing absorption (Asada+2024) |
| `cgm_A` | Float | `3.5918` | Sigmoid amplitude for HI column density evolution |
| `cgm_a` | Float | `1.8414` | Sigmoid slope for HI column density evolution |
| `cgm_c` | Float | `18.001` | Sigmoid baseline for log10(N\_HI) |

The CGM model adds damping wing absorption redward of Lyman-alpha at z >= 6, applied multiplicatively on top of the IGM attenuation. The HI column density evolves with redshift as:

```
log10(N_HI) = A / (1 + exp(-a * (z - 6))) + c
```

The default parameters are from Asada et al. (2024) and generally do not need to be changed.

## `[runtime]` -- Processing Options

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `chunked_processing` | Bool | `false` | Enable memory-efficient chunked processing |
| `target_memory_gb` | Float | `0.5` | Target memory usage per chunk in GB |
| `preserve_work_file` | Bool | `false` | Keep intermediate `.work.h5` file after completion |

See [Advanced Usage](@ref) for details on chunked processing and resume capability.

## `[translate]` -- Column Mapping

Maps filter names to flux/error columns in your input catalog.

```toml
[translate]
    f606w = {flux = 'f_f606w', error = 'e_f606w'}
    f814w = {flux = 'f_f814w', error = 'e_f814w'}
    f115w = {flux = 'f_f115w', error = 'e_f115w'}
```

The keys are filter names or nicknames (see [Filters](@ref) for the full list). The values must match column names in the input FITS catalog.

### Spectroscopic Redshift Column

When `use_zspec = true` in `[fitting]`, add a `zspec` entry:

```toml
[translate]
    zspec = 'z_spec'
```

This maps to the column containing spectroscopic redshifts. Objects with NaN, negative, or missing values are fit normally.
