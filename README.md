# Lazy.jl

Photometric redshift fitting in Julia. Largely based on [eazy-py](https://github.com/gbrammer/eazy-py). 
Designed to be fast, simple, and multithreaded. 

## Installation 

Lazy.jl is not a registered Julia package, so you need to install it from source.
The easiest way to do this is to use the included `install.sh` script. 
This will install the package and all dependencies, and set up a `lazy` 
executable in your path.

First, make sure you have `julia` installed. 
I recommend following the instructions [here](https://julialang.org/install/), 
using `juliaup`. Lazy.jl requires Julia 1.12 or later, and is tested on 1.12. 
If you installed with `juliaup`, you may need to manually install 1.12 and set 
it as your default julia channel using
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
This will install the `lazy` executable at `~/.julia/bin/`. You will need to add 
this to your path, by adding the following line to your `~/.bashrc` or 
`~/.zshrc` file:
```bash
export PATH="$PATH:$HOME/.julia/bin"
```
Then, run `source ~/.bashrc` or `source ~/.zshrc` to update your path.

Note that to update your installation (e.g. after pulling new changes from the 
repository), you can simply run `bash install.sh` again. This will update the 
package and all dependencies, and recompile the package if necessary.

## Usage

The `lazy` executable is a command line interface to the Lazy.jl package.

### Fitting

To fit a catalog, run `lazy fit`. The first time you run this, it will 
precompile the package (and all dependencies), which may take a few minutes. 
This is normal, and will not happen again unless you update the package or 
change the Julia version.
```
> lazy fit 
     __                        _ __
    / /  ____ _____ __  __    (_) /
   / /  / __ `/_  // / / /   / / /
  / /__/ /_/ / / // /_/ /   / / /
 /_____|__,_/ /___|__, (_)_/ /_/
                 /____/ /___/
v1.0.3 (1 threads)
====================================
Lazy.jl
usage: lazy fit -p <param_file> -t <nthreads>

  -p, --param     Path to the parameter file
  -t, --threads  Number of threads to use

```

### The parameter file

Lazy.jl uses a parameter file to specify all options related to fitting. 
The parameter file is a [TOML](https://toml.io/en/) file, which is a simple 
key-value format, divided into sections. 

Below is an example parameter file, which you can use as a starting point for 
your own. Run `lazy params > example_param.toml` to export the example param 
file to the current directory. 

```toml
################################################################################
# Example Lazy.jl parameter file, May 2025, v1.0.3
################################################################################

[io]
    # Path to the input catalog file, in fits format
    input_catalog = '/path/to/catalog/file.fits' 

    # Format for missing data in the input catalog (e.g. 'nan' or -99) 
    missing_data_format = 'nan' 

    # Name of the output fits file 
    output_file = '/path/to/lazy_results.fits'
    
    # Whether to output the redshift probability distribution for each object
    # If true, outputs to output_file in extension "PZ" 
    output_pz = true 
    
    # Whether to output the template set used, including IGM attenuation 
    # This is needed in order to reconstruct the best-fit SEDs later
    # If true, outputs to output_file in extension "TEMPL" 
    output_templates = true 

[fitting]
    # Template set to use for fitting. This should correspond to an entry in the
    # "template_directory.toml" file. Use `lazy list-templates` to see the 
    # available template sets. Alternatively, you can specify a custom template
    # set by providing a list of template files.
    
    template_set = 'sfhz' 
    # template_set = [
    #     'sfhz/corr_sfhz_13_bin3_av0.50.fits', 
    #     'sfhz/corr_sfhz_13_bin3_av0.01.fits', 
    #     ...
    # ]

    # Template error function to use (currently only 'template_error' supported)
    template_error = 'template_error' 

    # Scale factor for the template error function
    template_error_scale = 0.2 

    # IGM attenuation model to use (currently only 'inoue14' is supported)
    igm_model = 'inoue14'
    
    # Minimum number of valid (non-nan, snr>2) photometric points to run the fit
    nphot_min = 2 

    # Systematic error to add in quadrature to the photometric errors
    sys_err = 0.05

    # Redshift grid to use for fitting
    z_min = 0.0
    z_max = 20.0
    z_step = 0.01

[translate]
    # Tells the code which columns in the input catalog to use for the fluxes 
    # and errors in each band.The keys are the filter names, which correspond to
    # entries in the "filter_directory.toml" file. Use `lazy list-filters` to 
    # see the available filters. The entires should correspond to column names
    # in the input catalog

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

### Multithreading 

Lazy.jl takes advantage of multithreading in Julia, and will automatically use 
all available threads by default.
You can specify the number of threads to use with the `-t` or `--threads` option. 
Pass `--threads auto` to automatically use the maximum number of threads available. 

The key advantage of using multithreading in Julia (vs. multiprocessing in Python) 
is that it allows shared memory access. 


### Adding your templates, filters, etc. 

Lazy.jl takes in templates in the same format as `eazy` or `eazy-py`. To add templates, 
just copy them over to the `src/templates` directory in your Lazy.jl installation 
path. Then, add an entry to the `template_directory.toml` file, providing the paths
to the templates and a name for the template set. 

Similarly, to add filters, add the filter file to the `src/filter_files` directory, 
preferably in a subdirectory for the specific telescope/instrument. Filters 
should be provided as ascii files with wavelength in angstroms. Then, add an entry
to the `filter_directory.toml` file, specifying the path to the filter file and
any nicknames you want to use to reference that particular filter. For example, 
```toml
['spitzer/irac/ch2']
description = 'Spitzer/IRAC CH2 downloaded from SVO filter profile service'
nicknames = ['irac_ch2','IRAC_ch2','ch2','CH2','IRAC_45','irac_45']
wav_units = 'angstrom'
```

