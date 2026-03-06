# Output Formats

Lazy.jl can output results as FITS (default), HDF5, or both. The format is controlled by the `output_file` extension or the `output_format` parameter.

## FITS Output

The output is a multi-extension FITS file.

### SUMMARY Extension

The primary results table. Always present.

**Core columns:**

| Column | Type | Description |
|--------|------|-------------|
| `ID` | Int | Object identifiers from input catalog |
| `z_best` | Float | Best-fit redshift (minimum chi-squared) |
| `chi2` | Float | Best-fit chi-squared value |
| `z_l95` | Float | Lower 95% confidence bound (2.5th percentile of P(z)) |
| `z_l68` | Float | Lower 68% confidence bound (16th percentile) |
| `z_med` | Float | Median redshift (50th percentile of P(z)) |
| `z_u68` | Float | Upper 68% confidence bound (84th percentile) |
| `z_u95` | Float | Upper 95% confidence bound (97.5th percentile) |
| `coeffs` | Float array | Best-fit template coefficients (length = number of templates) |

**Per-band photometry:** One column per filter with the best-fit model flux (e.g., `f606w`, `f115w`).

**P(z) summary statistics:**

| Column | Type | Description |
|--------|------|-------------|
| `Sz` | Float | Center of the unit-width redshift bin containing the most P(z) |
| `Pz_cen` | Float | Integrated P(z) within 0.15(1+z) of z\_best |
| `Pz_zgtrzb2` | Float | Integrated P(z > z\_best - 2) |
| `pz_gt{Z}` | Float | P(z > Z) for each integer Z in the grid (e.g., `pz_gt3`, `pz_gt5`) |
| `Pz{n}` | Float | Integrated P(z) in integer-centered unit bins |

**Conditional columns** (depending on configuration):

| Column | Present when | Description |
|--------|-------------|-------------|
| `z_spec` | `use_zspec = true` | Input spectroscopic redshift |
| `M_UV`, `M_U`, `M_V`, `M_J` | `output_restframe_mags = true` | Rest-frame absolute magnitudes |
| `z_best_lowz` | `output_forced_lowz = true` | Best-fit redshift from forced low-z fit |
| `chi2_lowz` | `output_forced_lowz = true` | Chi-squared from forced low-z fit |
| `delta_chi2` | `output_forced_lowz = true` | chi2\_lowz - chi2\_best |
| `z_l95_lowz` ... `z_u95_lowz` | `output_forced_lowz = true` | Low-z P(z) quantiles |
| `coeffs_lowz` | `output_forced_lowz = true` | Low-z template coefficients |

### PZ Extension

Present when `output_pz = true`. Contains the full redshift probability distribution for each object.

| Column | Type | Description |
|--------|------|-------------|
| `ID` | Int | Object ID (-1 for the first row, which stores the redshift grid) |
| `Pz` | Float array | P(z) values (first row contains the redshift grid values) |

The first row is a sentinel with `ID = -1` and `Pz` = the redshift grid. Subsequent rows contain the P(z) for each object.

### TEMPL Extension

Present when `output_templates = true`. Contains the IGM-attenuated template SEDs at each redshift, enabling full SED reconstruction.

| Column | Type | Description |
|--------|------|-------------|
| `z` | Float | Redshift grid (-1 for the first row, which stores the wavelength grid) |
| `{template_name}` | Float array | Template flux at each redshift (first row contains the common wavelength grid in Angstroms) |

## Reconstructing Best-Fit SEDs

With both `output_templates = true` and the SUMMARY coefficients, you can reconstruct the best-fit SED for any object:

```python
from astropy.io import fits
import numpy as np

ID = 1
lazy_file = 'results.fits'

# Load results
lazy_summary = fits.getdata(lazy_file, ext=1)
i = np.where(lazy_summary['ID'] == ID)[0][0]
z_best = lazy_summary['z_best'][i]
coeffs = lazy_summary['coeffs'][i]

# Load templates
lazy_templ = fits.getdata(lazy_file, ext=3)
izbest = np.argmin(np.abs(z_best - lazy_templ['z']))
template_names = [n for n in lazy_templ.dtype.names if n != 'z']
templates = np.array([lazy_templ[t][izbest] for t in template_names])

# Reconstruct SED
wavelength = lazy_templ[template_names[0]][0] * (1 + z_best)  # Observed wavelength
fnu = np.dot(templates.T, coeffs)
```

## HDF5 Output

When the output file ends in `.h5` or `.hdf5`, or when `output_format = 'hdf5'`:

```
results.h5
├── metadata/           # Run parameters and configuration
├── results/            # Main fitting results (zbest, chi2, coeffs, quantiles, ...)
├── photometry/         # Best-fit model photometry per band
├── pz/                 # P(z) distributions (if output_pz = true)
│   ├── chi2grid         # Chi-squared at each redshift (nz x nobj)
│   └── pz               # Normalized P(z) (nz x nobj)
└── templates/          # Template SEDs (if output_templates = true)
```

HDF5 output uses compression and chunking for efficient storage of large datasets. It is also the format used internally for chunked processing work files.
