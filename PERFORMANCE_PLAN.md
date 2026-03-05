# Lazy.jl Performance Optimization Plan

Reference case: 1M objects, 20 bands, z=0-20 step 0.01 (nz=2001), ~13 templates.

## Item 1: Cache filter data (CRITICAL)

**Problem**: `get_filter()` (io.jl:95-119) re-parses `filter_directory.toml` and re-reads
filter files from disk on every call. Called ~520K times during template grid build.

**Fix**: Load all needed filters once into a Dict cache. Avoid repeated TOML parsing and disk I/O.

**Impact**: 10-100x faster template grid build.

**Status**: [x] Committed (5786ee2)

---

## Item 2: Reorder template grid to (nband, ntempl, nz) (CRITICAL)

**Problem**: `templgrid[ntempl, nz, nband]` forces non-contiguous copies when slicing
`templgrid[:,i,:]` in the inner fitting loop (2B allocations). Also requires manual transpose.

**Fix**: Reshape to `(nband, ntempl, nz)`. Slice `@view templgrid[:,:,i]` is contiguous and
already in the right orientation for NNLS. Update all code that reads/writes templgrid.

**Impact**: 2-3x faster fitting loop. Eliminates 2B allocations + transpose ops.

**Status**: [x] Committed (c20aaa3)

---

## Item 3: Hoist interpolation outside band loop (CRITICAL)

**Problem**: `linear_interpolation(wav_obs, templfnu_j)` is constructed once per
(template, redshift, band) = 520K times. Only depends on (template, redshift).

**Fix**: Build interpolation object once per (template, redshift), evaluate for each band.

**Impact**: 5-20x fewer interpolation constructions during grid build.

**Status**: [x] Committed (363b3d0)

---

## Item 4: In-place inner loop operations (HIGH)

**Problem**: `fit_single_object` creates ~8 temporary arrays per redshift step via
broadcasting (efnu_tot_j, valid, detect, snr_j). 2B iterations = ~16B temp arrays.

**Fix**: Pre-allocate all working arrays outside the redshift loop. Use `@.` for in-place
ops. Pre-compute `isfinite.(fnu_j)` once per object.

**Impact**: 1.5-2x faster fitting.

**Status**: [x] Committed (c7c0715)

---

## Item 5: NNLS workspace pre-allocation (HIGH)

**Problem**: `nonneg_lsq(templgrid_ij[valid,:], snr_j[valid])` copies into new arrays
every redshift step. The `[:]` forces another copy.

**Fix**: Pre-allocate workspace for valid-only NNLS problem. Copy valid rows into
pre-allocated buffers.

**Impact**: 1.3-1.5x faster fitting.

**Status**: [x] Implemented  [ ] Committed

---

## Item 6: P(z) quantiles via searchsorted (HIGH)

**Problem**: Each quantile creates a full `(nobj, nz)` temporary via
`abs.(cpz_chunk .- threshold)`. 5 quantiles = ~40 GB temporary allocations for 1M objects.

**Fix**: Use `searchsortedfirst` on the CDF per object. O(log n) per quantile, zero allocs.

**Impact**: 10x faster post-processing, eliminates ~40 GB temp allocations.

**Status**: [x] Implemented  [ ] Committed

---

## Item 7: Replace argmin(abs.()) with searchsorted (MEDIUM)

**Problem**: `argmin(abs.(zgrid .- value))` is O(n) with allocation. Used 9 times, some
in hot loops (lines 832, 961, 992, 1009 of fitting.jl).

**Fix**: Replace with `searchsortedlast` / `searchsortedfirst` (O(log n), zero alloc).

**Impact**: 1.1x overall.

**Status**: [ ] Planned  [ ] Implemented  [ ] Committed

---

## Item 8: chi2grid layout (nz, nobj) for cache efficiency (MEDIUM)

**Problem**: `chi2grid[nobj, nz]` means each object's P(z) row strides by nobj in memory.
For 1M objects, rows are 8 MB apart -> constant cache misses during P(z) computation.

**Fix**: Store as `(nz, nobj)` so each object's chi2 vector is contiguous. Update cumsum,
exp, eachrow -> eachcol, and HDF5 write patterns.

**Impact**: 2-5x faster P(z) computation.

**Status**: [x] Implemented  [ ] Committed

---

## Item 9: Pre-computed filter integration weights (MEDIUM)

**Problem**: Per (template, redshift, band), the code interpolates the template onto the
filter grid then integrates. Could pre-compute filter weights for dot-product integration.

**Fix**: For each filter, pre-compute normalized weights. Integration becomes a dot product.

**Impact**: 1.5-2x faster grid build.

**Status**: [ ] Planned  [ ] Implemented  [ ] Committed

---

## Item 10: Float32 template grid (LOW-MEDIUM)

**Problem**: Template grid uses Float64 but Float32 is sufficient for photometric precision.
Halves memory bandwidth.

**Fix**: Change `templgrid` to `Float32`. Adjust NNLS call and downstream code.

**Impact**: 1.1-1.2x faster fitting via better cache utilization.

**Status**: [ ] Planned  [ ] Implemented  [ ] Committed
