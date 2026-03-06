# Filters

Lazy.jl includes a comprehensive library of filter transmission curves organized by telescope and instrument. Use `lazy list-filters` to see all available filters.

## Using Filters

In the `[translate]` section, reference filters by their full path or any of their nicknames:

```toml
[translate]
    # These are all equivalent ways to reference HST/ACS F606W:
    f606w = {flux = 'f_f606w', error = 'e_f606w'}
    # F606W = {flux = 'f_f606w', error = 'e_f606w'}
    # hst_acs_f606w = {flux = 'f_f606w', error = 'e_f606w'}
    # hst/acs/f606w = {flux = 'f_f606w', error = 'e_f606w'}
```

Nickname matching is case-sensitive. Each filter has multiple nickname variants to accommodate different naming conventions.

## Available Filters

### JWST

**NIRCam Wide**

| Filter | Nicknames |
|--------|-----------|
| `jwst/nircam/f070w` | `f070w`, `F070W`, `nircam_f070w` |
| `jwst/nircam/f090w` | `f090w`, `F090W`, `nircam_f090w` |
| `jwst/nircam/f115w` | `f115w`, `F115W`, `nircam_f115w` |
| `jwst/nircam/f150w` | `f150w`, `F150W`, `nircam_f150w` |
| `jwst/nircam/f200w` | `f200w`, `F200W`, `nircam_f200w` |
| `jwst/nircam/f277w` | `f277w`, `F277W`, `nircam_f277w` |
| `jwst/nircam/f356w` | `f356w`, `F356W`, `nircam_f356w` |
| `jwst/nircam/f444w` | `f444w`, `F444W`, `nircam_f444w` |

**NIRCam Medium**

| Filter | Nicknames |
|--------|-----------|
| `jwst/nircam/f140m` | `f140m`, `F140M`, `nircam_f140m` |
| `jwst/nircam/f162m` | `f162m`, `F162M`, `nircam_f162m` |
| `jwst/nircam/f182m` | `f182m`, `F182M`, `nircam_f182m` |
| `jwst/nircam/f210m` | `f210m`, `F210M`, `nircam_f210m` |
| `jwst/nircam/f250m` | `f250m`, `F250M`, `nircam_f250m` |
| `jwst/nircam/f300m` | `f300m`, `F300M`, `nircam_f300m` |
| `jwst/nircam/f335m` | `f335m`, `F335M`, `nircam_f335m` |
| `jwst/nircam/f360m` | `f360m`, `F360M`, `nircam_f360m` |
| `jwst/nircam/f410m` | `f410m`, `F410M`, `nircam_f410m` |
| `jwst/nircam/f430m` | `f430m`, `F430M`, `nircam_f430m` |
| `jwst/nircam/f460m` | `f460m`, `F460M`, `nircam_f460m` |
| `jwst/nircam/f480m` | `f480m`, `F480M`, `nircam_f480m` |

**NIRCam Narrow**

| Filter | Nicknames |
|--------|-----------|
| `jwst/nircam/f164n` | `f164n`, `F164N`, `nircam_f164n` |
| `jwst/nircam/f187n` | `f187n`, `F187N`, `nircam_f187n` |
| `jwst/nircam/f212n` | `f212n`, `F212N`, `nircam_f212n` |
| `jwst/nircam/f323n` | `f323n`, `F323N`, `nircam_f323n` |
| `jwst/nircam/f405n` | `f405n`, `F405N`, `nircam_f405n` |
| `jwst/nircam/f466n` | `f466n`, `F466N`, `nircam_f466n` |
| `jwst/nircam/f470n` | `f470n`, `F470N`, `nircam_f470n` |

**MIRI**

| Filter | Nicknames |
|--------|-----------|
| `jwst/miri/f560w` | `f560w`, `F560W`, `miri_f560w` |
| `jwst/miri/f770w` | `f770w`, `F770W`, `miri_f770w` |
| `jwst/miri/f1000w` | `f1000w`, `F1000W`, `miri_f1000w` |
| `jwst/miri/f1130w` | `f1130w`, `F1130W`, `miri_f1130w` |
| `jwst/miri/f1280w` | `f1280w`, `F1280W`, `miri_f1280w` |
| `jwst/miri/f1500w` | `f1500w`, `F1500W`, `miri_f1500w` |
| `jwst/miri/f1800w` | `f1800w`, `F1800W`, `miri_f1800w` |
| `jwst/miri/f2100w` | `f2100w`, `F2100W`, `miri_f2100w` |
| `jwst/miri/f2550w` | `f2550w`, `F2550W`, `miri_f2550w` |

### HST

**ACS**

| Filter | Nicknames |
|--------|-----------|
| `hst/acs/f435w` | `f435w`, `F435W`, `hst_f435w`, `hst_acs_f435w` |
| `hst/acs/f606w` | `f606w`, `F606W`, `hst_f606w`, `hst_acs_f606w` |
| `hst/acs/f775w` | `f775w`, `F775W`, `hst_f775w`, `hst_acs_f775w` |
| `hst/acs/f814w` | `f814w`, `F814W`, `hst_f814w`, `hst_acs_f814w` |
| `hst/acs/f850lp` | `f850lp`, `F850LP`, `f850l`, `hst_f850lp`, `hst_acs_f850lp` |

**WFC3**

| Filter | Nicknames |
|--------|-----------|
| `hst/wfc3/f098m` | `f098m`, `F098M`, `hst_f098m`, `hst_wfc3_f098m` |
| `hst/wfc3/f105w` | `f105w`, `F105W`, `hst_f105w`, `hst_wfc3_f105w` |
| `hst/wfc3/f125w` | `f125w`, `F125W`, `hst_f125w`, `hst_wfc3_f125w` |
| `hst/wfc3/f140w` | `f140w`, `F140W`, `hst_f140w`, `hst_wfc3_f140w` |
| `hst/wfc3/f160w` | `f160w`, `F160W`, `hst_f160w`, `hst_wfc3_f160w` |

### Ground-Based

**Subaru HSC**

| Filter | Nicknames |
|--------|-----------|
| `subaru/hsc/g` | `hsc_g`, `HSC-g`, `HSC_g` |
| `subaru/hsc/r` | `hsc_r`, `HSC-r`, `HSC_r` |
| `subaru/hsc/i` | `hsc_i`, `HSC-i`, `HSC_i` |
| `subaru/hsc/z` | `hsc_z`, `HSC-z`, `HSC_z` |
| `subaru/hsc/y` | `hsc_y`, `HSC-y`, `HSC_y` |

**DECam**

| Filter | Nicknames |
|--------|-----------|
| `DECam/decam_u` | `decam_u`, `DECam_u` |
| `DECam/decam_g` | `decam_g`, `DECam_g` |
| `DECam/decam_r` | `decam_r`, `DECam_r` |
| `DECam/decam_i` | `decam_i`, `DECam_i` |
| `DECam/decam_z` | `decam_z`, `DECam_z` |
| `DECam/decam_y` | `decam_y`, `DECam_y` |

**CFHT**

| Filter | Nicknames |
|--------|-----------|
| `cfht/cfht_u` | `cfht_u`, `CFHT_u` |
| `cfht/cfht_wircam_j` | `cfht_j`, `CFHT_J`, `cfht_wircam_j` |
| `cfht/cfht_wircam_ks` | `cfht_ks`, `CFHT_Ks`, `cfht_wircam_ks` |

**UltraVISTA**

| Filter | Nicknames |
|--------|-----------|
| `uvista/Y` | `UVISTA-Y`, `uvista_Y`, `UltraVISTA_Y` |
| `uvista/J` | `UVISTA-J`, `uvista_J`, `UltraVISTA_J` |
| `uvista/H` | `UVISTA-H`, `uvista_H`, `UltraVISTA_H` |
| `uvista/Ks` | `UVISTA-Ks`, `uvista_Ks`, `UltraVISTA_Ks` |

### Space Telescopes

**Euclid**

| Filter | Nicknames |
|--------|-----------|
| `euclid/VIS` | `VIS`, `vis`, `euclid_vis` |
| `euclid/NISP_Y` | `NISP_Y`, `nisp_y`, `euclid_NISP_Y` |
| `euclid/NISP_J` | `NISP_J`, `nisp_j`, `euclid_NISP_J` |
| `euclid/NISP_H` | `NISP_H`, `nisp_h`, `euclid_NISP_H` |

**Spitzer IRAC**

| Filter | Nicknames |
|--------|-----------|
| `spitzer/irac/ch1` | `ch1`, `CH1`, `irac_ch1`, `IRAC_36` |
| `spitzer/irac/ch2` | `ch2`, `CH2`, `irac_ch2`, `IRAC_45` |

### Rest-Frame Filters

These are used internally for rest-frame absolute magnitude computation (see [Configuration](@ref)).

| Filter | Nickname | Description |
|--------|----------|-------------|
| `restframe/tophat_uv` | `rf_uv` | 1450-1550 A tophat for M\_UV |
| `restframe/bessell_u` | `rf_u` | Bessell U for M\_U |
| `restframe/bessell_v` | `rf_v` | Bessell V for M\_V |
| `restframe/twomass_j` | `rf_j` | 2MASS J for M\_J |

## Adding Custom Filters

1. Place your filter file in `src/filter_files/telescope/instrument/`:
   ```
   src/filter_files/my_telescope/my_filter.dat
   ```
   The file should be two-column ASCII: wavelength (Angstroms) and transmission.

2. Add an entry to `src/filter_files/filter_directory.toml`:
   ```toml
   ['my_telescope/my_filter']
   description = 'My custom filter'
   nicknames = ['my_filter', 'MY_FILTER']
   wav_units = 'angstrom'
   ```

3. Verify with `lazy list-filters`.

Contributions of new filters are very welcome! If you've added filter curves for an instrument or survey not yet included, please open a pull request so others can use them too. See [Development](@ref) for the contributing workflow.
