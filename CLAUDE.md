# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Lazy.jl is a Julia package for fast photometric redshift fitting, designed as a multithreaded alternative to eazy-py. It uses spectral energy distribution (SED) template matching with IGM attenuation models to estimate redshifts from photometric data.

## Key Commands

### Installation & Setup
- `bash install.sh` - Install the package and create `lazy` executable in `~/.julia/bin/`
- Requires Julia 1.12+

### Running the CLI
- `lazy fit -p <param_file>` - Run photometric redshift fitting
- `lazy fit -p <param_file> -t <nthreads>` - Run with specific thread count  
- `lazy list-templates` - Show available template sets
- `lazy list-filters` - Show available filter transmission curves
- `lazy params > example.toml` - Export example parameter file

### Testing & Development
- `julia --project=. -e "using Pkg; Pkg.test()"` - Run test suite
- `julia --project=.` - Enter Julia REPL with project environment

## Architecture

### Core Structure
- **Entry Point**: `src/Lazy.jl` - Module definition with CLI entry point (`julia_main()`)
- **Main Logic**: `src/main.jl` - Contains all photometric fitting algorithms and CLI parsing
- **Binary**: `bin/lazy` - Bash wrapper that calls Julia with proper threading/environment

### Key Components

**Template System**: 
- Templates stored in `src/templates/` with metadata in `template_directory.toml`
- Supports both ASCII (.dat) and FITS formats
- Can be redshift-dependent (2D flux arrays) or redshift-independent

**Filter System**:
- Filter transmission curves in `src/filter_files/` organized by telescope/instrument
- Metadata in `filter_directory.toml` with nickname support for easy referencing

**Fitting Algorithm** (`src/main.jl:fit()`):
1. Loads photometric catalog and builds template grid with IGM attenuation
2. Pre-computes template errors for all redshift/band combinations (memory optimization)
3. Performs by-object fitting using non-negative least squares (`nonneg_lsq`)
4. Outputs best-fit parameters, redshift probability distributions, and templates

**Data Flow**:
- Input: FITS catalog + TOML parameter file
- Processing: Template interpolation → IGM application → Photometric integration
- Output: Multi-extension FITS file (SUMMARY, PZ, TEMPL extensions)

### Memory Optimization Features
- Pre-computed template error grids to avoid repeated calculations
- By-object processing to minimize memory usage for large catalogs
- Elimination of large coefficient arrays (only store best-fit coefficients)

### Threading Strategy
- Uses Julia's native multithreading (`@threads`) for template grid construction
- Per-object fitting loop is threaded for scalability
- Thread count controlled via `-t` flag or `JULIA_NUM_THREADS`

## Configuration

### Parameter Files (TOML format)
Required sections:
- `[io]` - Input/output file paths and formats
- `[fitting]` - Template sets, redshift grids, error models
- `[translate]` - Maps filter names to catalog column names

### Adding Templates/Filters
1. Copy files to appropriate subdirectories in `src/templates/` or `src/filter_files/`
2. Update corresponding `.toml` directory file with metadata
3. Templates must follow eazy format; filters are ASCII wavelength/transmission

## Error Handling
- Custom `LazyError` exception type for user-facing errors
- Graceful handling of missing data and insufficient photometry
- Validation of parameter files and data integrity

## Python Integration
- Uses PyCall.jl to interface with astropy for FITS I/O via `src/writedata.py`
- Python dependencies automatically managed through Conda.jl