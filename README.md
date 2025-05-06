# Lazy.jl

[![Build Status](https://github.com/hollisakins/Lazy.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/hollisakins/Lazy.jl/actions/workflows/CI.yml?query=branch%3Amain)

Photometric redshift fitting in Julia. Largely based on [eazy-py](https://github.com/gbrammer/eazy-py). Designed to be fast, simple, and multithreaded. 

## Installation 

Lazy.jl is not a registered Julia package, so you need to install it from source.
The easiest way to do this is to use the included `install.sh` script. 
This will install the package and all dependencies, and set up a `lazy` executable in your path.

First, make sure you have `julia` installed. I recommend following the instructions [here](https://julialang.org/install/), using `juliaup`. 
Lazy.jl requires Julia 1.12 or later, and is tested on 1.12. If you installed with `juliaup`, you may need to manually install 1.12 and set it as your default julia channel using
```bash
juliaup add 1.12
juliaup default 1.12
```

Then, clone the repository and run the install script:
```
git clone https://github.com/hollisakins/Lazy.jl.git
cd Lazy.jl
bash install.sh
```
This will install the `lazy` executable at `~/.julia/bin/`. You will need to add this to your path, by adding the following line to your `~/.bashrc` or `~/.zshrc` file:
```bash
export PATH="$PATH:$HOME/.julia/bin"
```
Then, run `source ~/.bashrc` or `source ~/.zshrc` to update your path.

Note that to update your installation (e.g. after pulling new changes from the repository), you can simply run `bash install.sh` again. This will update the package and all dependencies, and recompile the package if necessary.

## Usage

The `lazy` executable is a command line interface to the Lazy.jl package.

### Fitting

To fit a catalog, run `lazy fit`. The first time you run this, it will precompile the package (and all dependencies), which may take a few minutes. This is normal, and will not happen again unless you update the package or change the Julia version.
```
> lazy fit 
     __                        _ __
    / /  ____ _____ __  __    (_) /
   / /  / __ `/_  // / / /   / / /
  / /__/ /_/ / / // /_/ /   / / /
 /_____|__,_/ /___|__, (_)_/ /_/
                 /____/ /___/
v1.0.0-DEV (1 threads)
====================================
Lazy.jl
usage: lazy fit -p <param_file> -t <nthreads>

  -p, --param     Path to the parameter file
  -t, --threads  Number of threads to use

```

### Multithreading 

Lazy.jl takes advantage of multithreading in Julia, and will automatically use all available threads by default.
You can specify the number of threads to use with the `-t` or `--threads` option. 
Pass `--thread auto` to automatically use the maximum number of threads available. 

The key advantage of using multithreading in Julia (vs. multiprocessing in Python) is that it allows shared memory access. 

### The parameter file

Lazy.jl uses a parameter file to specify all options related to fitting. 
The parameter file is a [TOML](https://toml.io/en/) file, which is a simple key-value format, dividedd into sections. 

Example parameter file:
```toml
[io]
    input_catalog = 'example_input_catalog.fits' # Path to the input FITS catalog
    intput_missing_data = 'nan' # Format for missing data in the input catalog (e.g. 'nan' or -99) (currently not implemented)

    output_file = 'example_lazy_output.fits' # Name of the output file (with extension)
    output_pz = true # Output the redshift probability distribution for each object
    output_templates = true # Output the best-fit template for each object

[fitting]
    # Specify the template set to use for fitting
    # The template set must be in the Lazy.jl template list (see `lazy list-templates`)
    template_set = 'fsps_full' 

    # Template error function 
    template_error = 'template_error' # template error function to use
    template_error_scale = 0.2 # scale factor for the template error function

    # IGM model
    igm = 'igm_inoue14_grid' 
    
    # Redshift grid
    z_min = 0.0 
    z_max = 15.0
    z_step = 0.01

    # Fitting options
    sys_err = 0.05 # systematic error (% of flux) to add to the flux errors
    nphot_min = 2 # minimum number of valid photometric points to fit an object

[translate]
    # Specify the columns in the input catalog that correspond to the flux and error for each filter
    # The keys are filter names, which must be in the Lazy.jl filter list (see `lazy list-filters`)
    f115w = {flux = 'f_115', error = 'e_115'}
    f150w = {flux = 'f_150', error = 'e_150'}
    f200w = {flux = 'f_200', error = 'e_200'}
    f277w = {flux = 'f_277', error = 'e_277'}
    f356w = {flux = 'f_356', error = 'e_356'}
    f444w = {flux = 'f_444', error = 'e_444'}
```


