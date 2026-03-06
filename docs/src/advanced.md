# Advanced Usage

## Chunked Processing

For large catalogs (>100k objects) or fine redshift grids, enable chunked processing to keep memory usage bounded:

```toml
[runtime]
    chunked_processing = true
    target_memory_gb = 1.0
```

Lazy.jl automatically calculates the optimal chunk size based on your target memory and the dimensions of the problem (number of redshifts, bands, and templates). Multithreading operates within each chunk. Results are streamed to an HDF5 work file and converted to the final output format on completion.

## Resume Capability

If a job is interrupted (Ctrl-C, crash, etc.), Lazy.jl can resume from where it left off. Progress is saved in a `.work.h5` file alongside the output.

```bash
# If the previous run was interrupted:
lazy fit -p params.toml
# Will prompt: "Found incomplete run (75% complete), resume? [Y/n]"
```

For automated pipelines, use the `--yes` flag to skip the prompt and automatically resume:

```bash
lazy fit -p params.toml -y
```

To discard the work file and start fresh, delete the `.work.h5` file or choose "restart" when prompted.

The `preserve_work_file` option controls whether the work file is kept after successful completion:

```toml
[runtime]
    preserve_work_file = false   # Delete work file after completion (default)
```

## Template Grid Caching

Template grids (the pre-computed template photometry at every redshift) are cached to disk by default. This avoids rebuilding the grid when re-running with the same parameters.

- **First run**: Builds and caches the grid (~1 minute for a typical setup)
- **Subsequent runs**: Loads from cache (~2 seconds)

The cache key includes the template set, redshift grid, filter set, IGM model, template error parameters, and CGM settings. Any change to these parameters triggers a rebuild.

To manage the cache:

```bash
lazy cache-clear    # Remove all cached grids
```

Disable caching entirely with:

```toml
[fitting]
    template_cache = false
```

## Spectroscopic Redshifts

When spectroscopic redshifts are available for some objects, Lazy.jl can fix the fit at those redshifts:

```toml
[fitting]
    use_zspec = true

[translate]
    zspec = 'z_spec'
```

Objects with a valid z\_spec (positive, finite) are fit only at the nearest redshift grid point. Objects without z\_spec (NaN, negative, or missing column) are fit normally across the full grid. The output includes a `z_spec` column.

## Rest-Frame Absolute Magnitudes

Lazy.jl can compute rest-frame absolute magnitudes using the best-fit SED:

```toml
[io]
    output_restframe_mags = true
    flux_units = 'uJy'
    H0 = 70.0       # optional, default 70.0
    Om = 0.3         # optional, default 0.3
```

This outputs four additional columns: `M_UV` (1500 A tophat), `M_U` (Bessell U), `M_V` (Bessell V), and `M_J` (2MASS J). The computation uses a flat LCDM cosmology for the distance modulus and applies a K-correction of 2.5 log10(1+z).

Supported flux units and their AB zero points:

| Unit | Zero Point |
|------|-----------|
| `uJy` | 23.9 |
| `nJy` | 31.4 |
| `Jy` | 8.9 |
| `cgs` | -48.6 |

## Forced Low-Redshift Fitting

For high-redshift candidate validation, Lazy.jl can perform a parallel low-z fit:

```toml
[io]
    output_forced_lowz = true
    forced_lowz_zmax = 7.0
```

This runs a second fit for each object restricted to z < `forced_lowz_zmax`. The output includes:

- `z_best_lowz`, `chi2_lowz`: Best-fit low-z results
- `delta_chi2`: chi2\_lowz - chi2\_best (larger values indicate stronger evidence for the high-z solution)
- Low-z P(z) quantiles and template coefficients
- Low-z model photometry for each band

## CGM Damping Wing

The CGM Lyman-alpha damping wing model (Asada et al. 2024) adds absorption redward of Lyman-alpha at z >= 6. This is enabled by default:

```toml
[fitting]
    add_cgm = true
```

The model parameterizes the HI column density evolution as a sigmoid:

```
log10(N_HI) = A / (1 + exp(-a * (z - 6))) + c
```

The default parameters (`cgm_A = 3.5918`, `cgm_a = 1.8414`, `cgm_c = 18.001`) are from Asada et al. (2024) and generally do not need to be changed. The absorption is applied multiplicatively on top of the Inoue+2014 IGM attenuation.

## Performance and Scaling

### Threading

Scaling with thread count (1M objects, 10 bands, 20 templates):

| Threads | Time | Speedup | Efficiency |
|---------|------|---------|------------|
| 1 | 45 min | 1.0x | 100% |
| 4 | 12 min | 3.8x | 95% |
| 8 | 6.5 min | 6.9x | 86% |
| 16 | 4.2 min | 10.7x | 67% |

### Memory Usage

| Dataset | Objects | Redshifts | In-Memory | Chunked |
|---------|---------|-----------|-----------|---------|
| Small | 10k | 200 | 0.02 GB | N/A |
| Medium | 100k | 600 | 0.45 GB | N/A |
| Large | 500k | 2000 | 3.7 GB | 1.0 GB |
| Very Large | 1M | 2000 | 7.5 GB | 1.0 GB |
| Extreme | 10M | 2000 | 75 GB | 1.0 GB |

### Optimization Tips

1. **Redshift step size**: A fine grid (0.001) uses ~100x more memory than a coarse grid (0.01). Start coarse and refine as needed.
2. **Chunked processing**: Enable for datasets with estimated memory >1 GB.
3. **Template caching**: Keep enabled (default) for significant speedups on repeated runs.
4. **Thread count**: Optimal is typically 0.5-1x your CPU core count. Diminishing returns above ~8 threads.
5. **Memory target**: Set `target_memory_gb` to match your available system RAM, leaving headroom for the OS.
