# Templates

Lazy.jl fits photometric data by finding the best linear combination of spectral energy distribution (SED) templates at each redshift. The choice of template set determines what types of galaxies the fitter can model.

## Built-in Template Sets

Use `lazy list-templates` to see available sets. Select a set in your parameter file with:

```toml
[fitting]
    template_set = 'sfhz'
```

### Base Sets

| Template Set | Templates | Attribution | Description |
|-------------|-----------|-------------|-------------|
| `fsps_full` | 12 | eazy, FSPS | Standard FSPS templates from eazy |
| `sfhz` | 13 | eazy, FSPS, Gabe Brammer | Star-formation-history templates with redshift-dependent SEDs |
| `jades` | 16 | Hainline+, Brammer+08, Coe+06, Erb+10 | JADES survey template set |

### Extended Sets

The base sets can be extended with additional template components:

- **Larson22**: High-redshift galaxy templates from Larson et al. (2022), designed for z > 6 candidates with young stellar populations and nebular emission
- **Schrodinger** (`_schr`): Placeholder template for ambiguous sources
- **LRD** (`_lrd`): Little Red Dot template for compact red sources at high redshift

| Template Set | Templates | Components |
|-------------|-----------|------------|
| `fsps_Larson22` | 15 | FSPS + Larson22 |
| `fsps_Larson22_schr` | 16 | FSPS + Larson22 + Schrodinger |
| `fsps_Larson22_lrd` | 16 | FSPS + Larson22 + LRD |
| `fsps_Larson22_schr_lrd` | 17 | FSPS + Larson22 + Schrodinger + LRD |
| `sfhz_Larson22` | 19 | SFHZ + Larson22 |
| `sfhz_Larson22_schr` | 20 | SFHZ + Larson22 + Schrodinger |
| `sfhz_schr` | 14 | SFHZ + Schrodinger |
| `sfhz_lrd` | 14 | SFHZ + LRD |
| `sfhz_schr_lrd` | 15 | SFHZ + Schrodinger + LRD |
| `sfhz_Larson22_lrd` | 19 | SFHZ + Larson22 + LRD |
| `sfhz_Larson22_schr_lrd` | 20 | SFHZ + Larson22 + Schrodinger + LRD |
| `jades_Larson22` | 19 | JADES + Larson22 |

## Custom Template Sets

Instead of a named set, you can provide a list of template file paths:

```toml
[fitting]
    template_set = [
        'sfhz/corr_sfhz_13_bin3_av0.50.fits',
        'sfhz/corr_sfhz_13_bin3_av0.01.fits',
        'custom/my_template.dat',
    ]
```

Paths are relative to the `src/templates/` directory.

## Template File Formats

Lazy.jl supports two template formats:

**ASCII (`.dat`, `.sed`)**: Two-column files with wavelength (Angstroms) and flux. These are redshift-independent -- the same SED is used at all redshifts.

**FITS**: Contains a 2D flux array that varies with redshift. These are redshift-dependent templates where the SED shape changes across the fitting grid. The SFHZ templates use this format.

## Adding Custom Templates

1. Place your template files in a subdirectory of `src/templates/`:
   ```
   src/templates/my_set/template1.dat
   src/templates/my_set/template2.dat
   ```

2. Add an entry to `src/templates/template_directory.toml`:
   ```toml
   [my_set]
       files = [
           'my_set/template1.dat',
           'my_set/template2.dat',
       ]
       attribution = "Your citation"
   ```

3. Verify with `lazy list-templates`.

Contributions of new template sets are very welcome! If you've added templates that may be useful to others, please open a pull request to include them in the package. See [Development](@ref) for the contributing workflow.
