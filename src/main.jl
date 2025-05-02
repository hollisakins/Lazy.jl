using TOML
using FITSIO
using NonNegLeastSquares, Glob
using DelimitedFiles
using ProgressBars
using Interpolations
using Trapz
using LinearAlgebra
using Base.Threads
using HDF5
using Comonicon
version = "0.1.0"
filterpath = @__DIR__() * "/filter_files/"
igmpath = @__DIR__() * "/igm_data/"
templatepath = @__DIR__() * "/templates/"

"""
Fit photometric redshifts with Lazy.jl

# Args

- `param_file`: path to TOML parameter file

"""


function panic(
    msg::String, err::Union{Exception, Nothing} = nothing,
    bt::Union{Vector{Base.StackFrame}, Nothing} = nothing
)
    printstyled(stderr, "ERROR: "; color = :red, bold = true)
    print(stderr, msg)
    if err !== nothing
        print(stderr, sprint_showerror(err))
    end
    if bt !== nothing
        Base.show_backtrace(stderr, bt)
    end
    println(stderr)
    global errno = 1
    return errno
end

function print_help(io)
    printstyled(io, "Lazy.jl \n", bold = true)
    println("Usage: lazy -p <param_file> [-o <output_file>]")
    println("Options:")
    println("  -p, --param   Path to the parameter file")
    println("  -o, --output  Path to the output file")
end

function print_ascii(io)
    println("     __                        _ __ ")
    println("    / /  ____ _____ __  __    (_) / ")
    println("   / /  / __ `/_  // / / /   / / /  ")
    println("  / /__/ /_/ / / // /_/ /   / / /   ")
    println(" /_____|__,_/ /___|__, (_)_/ /_/    ")
    println("                 /____/ /___/       ")
    println("v$version")
    println("====================================")
end

function main(argv)::Cint
    
    io = stdout
    print_ascii(io)

    if argv == []
        print_help(io)
        return 1
    end

    # Parse command line arguments
    param = nothing
    while length(argv) > 0
        x = popfirst!(argv)
        if x == "-h" || x == "--help"
            print_help(io)
            return 1
        elseif x == "-p" || x == "--param"
            if length(argv) < 1
                return panic("expected parameter file argument after `-p`")
            end
            param = popfirst!(argv)
        else
            # Argument not recognized
            return panic("unrecognized argument: $x")
        end
    end
    if param == nothing
        return panic("parameter file not specified. Use -p <param_file>")
    end

    # Load in TOML parameter file
    if !isfile(param)
        return panic("parameter file $param not found.")
    end
    if !endswith(param, ".toml")
        return panic("parameter file $param must have a .toml extension.")
    end
    println("Loading parameter file: $param")
    param = TOML.parsefile(param)
    

    # Attempt to read in the input catalog file
    if haskey(param, "io")
        io = param["io"]
        if haskey(io, "input_catalog")
            input_catalog = io["input_catalog"]
            println("Input catalog: $input_catalog")
            # Read the input catalog
            cat = FITS(input_catalog)
        else
            panic("parameter `input_catalog` not found in the parameter file.")
        end
    else
        panic("section `io` not found in the parameter file.")
    end

    IDs = read(cat[2], "ID")
    nobj = length(IDs)
    println("nobj = $nobj")
        
    println("====================================")
    if !haskey(param, "translate")
        error("section `translate` not found in the parameter file.")
    end
    
    translate = param["translate"]
    bands = sort(collect(keys(param["translate"])))
    nband = length(bands)
    println("nband = $nband")
    println("bands = $bands")

    # Check that all bands exist as valid filter files
    for band in bands
        filterfile = filterpath * band * ".dat"
        if !isfile(filterfile)
            panic("Filter file $filterfile does not exist.")
        end
    end

    # Load in the data
    fnu, efnu = load_data(cat, bands, translate)
    println("fnu = " * summary(fnu))
    println("efnu = " * summary(efnu))
    println("====================================")
    
    if !haskey(param, "fitting")
        panic("section `fitting` not found in the parameter file.")
    end

    fitting = param["fitting"]

    if !haskey(fitting, "nphot_min")
        panic("parameter `nphot_min` not found in the parameter file.")
    end
    nphot_min = fitting["nphot_min"]

    if !haskey(fitting, "sys_err")
        panic("parameter `sys_err` not found in the parameter file.")
    end
    sys_err = fitting["sys_err"]
    fnu, efnu = set_sys_err(fnu, efnu, sys_err)

    if !haskey(fitting, "z_min")
        panic("parameter `z_min` not found in the parameter file.")
    end
    if !haskey(fitting, "z_max")
        panic("parameter `z_max` not found in the parameter file.")
    end
    if !haskey(fitting, "z_step")
        panic("parameter `z_step` not found in the parameter file.")
    end

    z_min = fitting["z_min"]
    z_max = fitting["z_max"]
    z_step = fitting["z_step"]
    println("z_min = $z_min")
    println("z_max = $z_max")
    println("z_step = $z_step")

    zgrid = collect(range(z_min, stop=z_max, step=z_step))
    nz = length(zgrid)

    println("====================================")

    if !haskey(fitting, "template_set")
        error("parameter `template_set` not found in the parameter file.")
    end

    template_set = fitting["template_set"]
    println("Template set: $template_set")
 
    templates = glob("*", templatepath * template_set)
    ntempl = length(templates)
    s = 0
    for (i, templ) in enumerate(templates)
        t = readdlm(templ)
        wave = t[:,1]
        flux = t[:,2]
        if i == 1
            s = size(t)[1]
        else
            if size(t)[1] != s
                error("Template $i: $templ has different number of rows than the first template.")
            end
        end
        println("Template $i: $templ $s")
        
    end
    
    templwav = readdlm(templates[1])[:,1]
    templfnu = zeros(s, ntempl)
    for (i, templ) in enumerate(templates)
        t = readdlm(templ)
        wave = t[:,1]
        flux = t[:,2]
        templfnu[:, i] =  flux .* wave .^ 2
    end

    templgrid = zeros(ntempl, nz, nband)
    println("Building template grid: " * summary(templgrid))
    
    if !haskey(fitting, "igm")
        error("parameter `igm` not found in the parameter file.")
    end
    
    h5open(igmpath * fitting["igm"] * ".hdf5", "r") do igm_file
        igm_redshifts = igm_file["redshifts"][:]
        igm_wavelengths = igm_file["wavelengths"][:]
        igm_transmission = igm_file["transmission"][:,:]
        

        iter = ProgressBar(1:ntempl)
        for i in iter
            @threads for j in 1:nz
                z = zgrid[j]

                templfnu_i = templfnu[:, i]
                wav_obs = templwav .* (1+z)
                
                # Interpolate the IGM transmission at this redshift
                # first interpolate over the redshift
                iz_up = searchsortedfirst(igm_redshifts, z)
                transmission = igm_transmission[iz_up,:]

                idx  = searchsortedfirst(templwav, 1215.67)
                idx_igm = searchsortedfirst(igm_wavelengths, 1215.67)
                interp = linear_interpolation([0.0;igm_wavelengths[1:idx_igm-1]], [0.0;transmission[1:idx_igm-1]], extrapolation_bc=Flat())
                y1 = interp(templwav[1:idx-1])
                interp = linear_interpolation([igm_wavelengths[idx_igm:end];1226.0], [transmission[idx_igm:end];1.0], extrapolation_bc=Flat())
                y2 = interp(templwav[idx:end])
                transmission  = [y1; y2]

                templfnu_i = templfnu_i .* transmission
                
                for k in 1:nband
                    
                    band = bands[k]
                    filt = readdlm(filterpath * band * ".dat")
                    fwav = filt[:,1]
                    ftrans = filt[:,2]
                    nu = 1 ./ fwav
                    interp = linear_interpolation(wav_obs, templfnu_i)
                    fnu_interp = interp(fwav)
                
                    result = trapz(nu, fnu_interp .* ftrans ./ nu) / trapz(nu, ftrans ./ nu)
                    templgrid[i,j,k] = result
                end
            end
        end
    end
   
    if !haskey(fitting, "template_error")
        error("parameter `template_error` not found in the parameter file.")
    end
   
    if !haskey(fitting, "template_error_scale")
        error("parameter `template_error_scale` not found in the parameter file.")
    end

    template_error = fitting["template_error"]
    template_error_scale = fitting["template_error_scale"]
    println("Template error: $template_error")

    tef = readdlm(template_error * ".dat")
    tef_x = tef[:,1]
    tef_y = tef[:,2]
    tef_xmin = minimum(tef_x[tef_y .> 0])
    tef_xmax = maximum(tef_x[tef_y .> 0])
    tef_clip_min = tef_y[tef_y .> 0][1]
    tef_clip_max = tef_y[tef_y .> 0][end]
    println("TEF x range: $tef_xmin - $tef_xmax")
    println("TEF y range: $tef_clip_min - $tef_clip_max")

    

    # Compute pivot wavelengths
    pivot_wavs = zeros(nband)
    for k in 1:nband
        band = bands[k]
        filt = readdlm(filterpath * band * ".dat")
        fwav = filt[:,1]
        ftrans = filt[:,2]
        pivot_wavs[k] = sqrt(trapz(fwav, ftrans) / trapz(fwav, ftrans ./ (fwav .^ 2)))
    end

    println("===========================")
    println("Fitting by redshift ")
    
    chi2grid = zeros(nobj,nz)
    
    iter = ProgressBar(1:nz)
    for i in iter
        
        templgrid_i = templgrid[:,i,:]
        chi2 = zeros(nobj)

        pivot_wavs_rest = pivot_wavs ./ (1+zgrid[i])
        tef_interp = linear_interpolation(tef_x, tef_y, extrapolation_bc=Flat())
        tefz = tef_interp(pivot_wavs_rest)
        tefz[pivot_wavs_rest .< tef_xmin] .= tef_clip_min
        tefz[pivot_wavs_rest .> tef_xmax] .= tef_clip_max
        tefz = tefz * template_error_scale

        @threads for j in 1:nobj
            fnu_j = fnu[j,:]
            efnu_j = efnu[j,:]
            efnu_tot_j = sqrt.( efnu_j .^ 2 + (tefz .* max.(fnu_j, 0.0)) .^ 2 )
            snr_j = fnu_j ./ efnu_tot_j
            
            #good = ~np.isnan(fnu_i) & ~np.isnan(efnu_tot_i)

            templgrid_ij = transpose(templgrid_i) ./ efnu_tot_j
            
            result = nonneg_lsq(templgrid_ij, snr_j ; alg=:nnls)[:]
            # println("Result: " * summary(result))

            fnu_mod_j = sum(result .* templgrid_i, dims=1)[:]
            chi2[j] = sum(((fnu_j .- fnu_mod_j) .^ 2) ./ efnu_tot_j .^ 2)
            
        end

        chi2grid[:,i] = chi2

    end
    
    zbest = zgrid[map(argmin, eachrow(chi2grid))]


    snr = fnu ./ efnu
    nphot = sum(snr .>= 2.0, dims=2)
    zbest[nphot .< nphot_min] .= -1

    

    return 0
end


function load_data(cat, bands, translate)
    """
    Load the photometric data from the catalog.
    """
    IDs = read(cat[2], "ID")
    nobj = length(IDs)
    nband = length(bands)
    fnu = zeros(nobj, nband)
    efnu = zeros(nobj, nband)
    for (i, band) in enumerate(bands)
        flux_col = translate[band]["flux"]
        err_col = translate[band]["error"]
        fnu_i = read(cat[2], flux_col)
        efnu_i = read(cat[2], err_col)
        fnu[:, i] = fnu_i
        efnu[:, i] = efnu_i
    end
    return fnu, efnu
end


function set_sys_err(fnu, efnu, sys_err)
    """
    Set the systematic error on the fluxes.
    """
    if sys_err == 0.0
        return fnu, efnu
    end
    nbands = size(fnu, 2)
    for i in 1:nbands
        fnu_i = fnu[:, i]
        efnu_i = efnu[:, i]
        efnu_i = sqrt.( efnu_i .^ 2 + (sys_err .* max.(fnu_i, 0.0)) .^ 2 )
        efnu[:, i] = efnu_i
    end
    return fnu, efnu
end

function list_filters()
    """
    List the available filters in the filter_files directory.
    """
    filter_files = glob("*.dat", filterpath)
    println("Available filters:")
    for file in filter_files
        println(file)
    end
end


function list_templates()
    """
    List the available templates in the template_set directory.
    """
    
    template_set = fitting["template_set"]
    println("Template set: $template_set")

    templates = glob(fitting["template_set"] * "/*")
    println("Available templates:")
    for templ in templates
        println(templ)
    end
end


Base.@main
