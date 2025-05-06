errno::Cint = 0
supports_color(io) = get(io, :color, false)

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
using PyCall



igmpath = @__DIR__() * "/igm_data/"
templatepath = @__DIR__() * "/templates/"
template_directory = TOML.parsefile(templatepath * "template_directory.toml")
template_set_names = sort(collect(keys(template_directory)))

filterpath = @__DIR__() * "/filter_files/"
filter_directory = TOML.parsefile(filterpath * "filter_directory.toml")
filter_nicknames = Dict{String, String}()
filter_names = sort(collect(keys(filter_directory)))

for key in filter_names
    for n in filter_directory[key]["nicknames"]
        filter_nicknames[n] = key
    end
end

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

function print_help(io, cmd::String = "main")
    printstyled(io, "Lazy.jl \n", bold = true)
    if cmd == "main"
        println("usage: lazy <command> <options>")
        println("")
        println("  fit <options>       Fit the data")
        println("  list-templates      List the available templates")
        println("  list-filters        List the available filters")
    elseif cmd == "fit"
        println("")
        println("  -p, --param     Path to the parameter file")
        println("  -t, -- threads  Number of threads to use")
    elseif cmd == "list-templates"
        println("usage: lazy list-templates")
    end

end

function print_ascii(io)
    nthreads = Threads.nthreads()
    println("     __                        _ __ ")
    println("    / /  ____ _____ __  __    (_) / ")
    println("   / /  / __ `/_  // / / /   / / /  ")
    println("  / /__/ /_/ / / // /_/ /   / / /   ")
    println(" /_____|__,_/ /___|__, (_)_/ /_/    ")
    println("                 /____/ /___/       ")
    if nthreads > 1
        println("v$version ($nthreads threads)")
    else
        println("v$version (1 thread)")
    end
    println("====================================")
end

function main(argv)
    check_version()

    # Reset errno
    global errno = 0
    
    io = stdout
    print_ascii(io)

    if argv == []
        print_help(io)
        return errno
    end

    # Parse command line arguments
    while length(argv) > 0
        x = popfirst!(argv)
        
        # Return the help message and exit
        if x == "-h" || x == "--help"
            print_help(io)
            return errno
        
        # 
        elseif x == "fit"
            if argv == []
                print_help(io, "fit")
                return errno
            end
        
            y = popfirst!(argv)

            if y == "-p" || y == "--param"
                if length(argv) < 1
                    return panic("expected parameter file argument after `-p`")
                end
                param = popfirst!(argv)
                return fit(param)
            elseif y == "-h" || y == "--help"
                print_help(io, "fit")
                return errno
            else
                return panic("unrecognized argument: $y")
            end


        elseif x == "list-templates"
            # List the available templates
            return list_templates()

        elseif x == "list-filters"
            # List the available templates
            return list_filters()
            
        else
            # Argument not recognized
            return panic("unrecognized argument: $x")
        end
    end
    if param == nothing
        return panic("parameter file not specified. Use -p <param_file>")
    end

end

function fit(param)

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
            panic("parameter `io.input_catalog` not found in the parameter file.")
        end

        if haskey(io, "output_file")
            output_file = io["output_file"]
        else
            panic("parameter `io.output_file` not found in the parameter file.")
        end

        if haskey(io, "output_pz")
            output_pz = io["output_pz"]
        else
            output_pz = false
        end
        if haskey(io, "output_templates")
            output_templ = io["output_templates"]
        else
            output_templ = false
        end
        output_file_pz = replace(output_file, ".fits" => "_pz.fits")
        output_file_templ = replace(output_file, ".fits" => "_templ.fits")
    else
        panic("section `io` not found in the parameter file.")
    end

    IDs = Int.(read(cat[2], "ID"))
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

    # Load in the data
    fnu, efnu = load_data(cat, bands, translate)
    println("fnu = " * summary(fnu))
    println("efnu = " * summary(efnu))

    if !haskey(io, "missing_data_format")
        panic("parameter `io.missing_data_format` not found in the parameter file.")
    end
    missing_data_format = io["missing_data_format"]
    if missing_data_format == "nan" || missing_data_format == "NaN"
        begin
        end
    elseif missing_data_format == "zero"
        fnu[efnu .= 0.0] .= NaN
        efnu[efnu .= 0.0] .= NaN
    elseif missing_data_format isa Number
        condition = (fnu .== missing_data_format) .| (efnu .== missing_data_format)
        fnu[condition] .= NaN
        efnu[condition] .= NaN
    else
        println("Warning: unrecognized missing data format: $missing_data_format, using input catalog as is.")
    end

        


    println("====================================")

    
    if !haskey(param, "fitting")
        panic("section `fitting` not found in the parameter file.")
    end

    fitting = param["fitting"]

    if !haskey(fitting, "nphot_min")
        panic("parameter `fitting.nphot_min` not found in the parameter file.")
    end
    nphot_min = fitting["nphot_min"]

    if !haskey(fitting, "sys_err")
        panic("parameter `fitting.sys_err` not found in the parameter file.")
    end
    sys_err = fitting["sys_err"]
    fnu, efnu = set_sys_err(fnu, efnu, sys_err)

    if !haskey(fitting, "z_min")
        panic("parameter `fitting.z_min` not found in the parameter file.")
    end
    if !haskey(fitting, "z_max")
        panic("parameter `fitting.z_max` not found in the parameter file.")
    end
    if !haskey(fitting, "z_step")
        panic("parameter `fitting.z_step` not found in the parameter file.")
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
        panic("parameter `template_set` not found in the parameter file.")
    end

    template_set = fitting["template_set"]
    if !haskey(template_directory, template_set)
        panic("Template set $template_set not found in the template directory.")
    end
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
        templ_shortname = basename(templ)
        println("Template $i: $templ_shortname ($s,)")
        
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
    
    igm_file = h5open(igmpath * fitting["igm"] * ".hdf5", "r")
    igm_redshifts = igm_file["redshifts"][:]
    igm_wavelengths = igm_file["wavelengths"][:]
    igm_transmission = igm_file["transmission"][:,:]
    close(igm_file)
    idx  = searchsortedfirst(templwav, 1215.67)
    idx_igm = searchsortedfirst(igm_wavelengths, 1215.67)

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

            interp = linear_interpolation([0.0;igm_wavelengths[1:idx_igm-1]], [0.0;transmission[1:idx_igm-1]], extrapolation_bc=Flat())
            y1 = interp(templwav[1:idx-1])
            interp = linear_interpolation([igm_wavelengths[idx_igm:end];1226.0], [transmission[idx_igm:end];1.0], extrapolation_bc=Flat())
            y2 = interp(templwav[idx:end])
            transmission  = [y1; y2]

            templfnu_i = templfnu_i .* transmission
            
            @threads for k in 1:nband
                
                band = bands[k]
                fwav, ftrans = get_filter(band)
                nu = 1 ./ fwav
                interp = linear_interpolation(wav_obs, templfnu_i)
                fnu_interp = interp(fwav)
            
                result = trapz(nu, fnu_interp .* ftrans ./ nu) / trapz(nu, ftrans ./ nu)
                templgrid[i,j,k] = result
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

    tef = readdlm(templatepath * "template_error/" * template_error * ".dat")
    tef_x = tef[:,1]
    tef_y = tef[:,2]
    tef_xmin = minimum(tef_x[tef_y .> 0])
    tef_xmax = maximum(tef_x[tef_y .> 0])
    tef_clip_min = tef_y[tef_y .> 0][1]
    tef_clip_max = tef_y[tef_y .> 0][end]
    println("TEF x range: $tef_xmin - $tef_xmax")
    println("TEF y range: $tef_clip_min - $tef_clip_max")
    pivot_wavs = get_pivot_wavelengths(bands)

    

    println("===========================")
    println("Fitting by redshift ")
    
    chi2grid = zeros(nobj,nz)
    coeffs = zeros(nobj,nz,ntempl)
    
    iter = ProgressBar(1:nz)
    for i in iter
        
        templgrid_i = templgrid[:,i,:]

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

            valid = isfinite.(fnu_j) .&  isfinite.(efnu_tot_j)
            detect = valid .& (snr_j .> 2.0)
            if sum(detect) < nphot_min
                #println("Object $j: not enough detections (nphot = $(sum(detect)))")
                chi2grid[j,i] = -1
                continue
            end

            templgrid_ij = transpose(templgrid_i) ./ efnu_tot_j
            
            result = nonneg_lsq(templgrid_ij[valid,:], snr_j[valid] ; alg=:nnls)[:]
            coeffs[j,i,:] = result

            # fnu_mod_j = sum(result .* templgrid_i, dims=1)[:]
            fnu_mod_j = transpose(templgrid_i) * result # order is important here
            chi2grid[j,i] = sum(((fnu_j .- fnu_mod_j) .^ 2) ./ efnu_tot_j .^ 2)
            
        end

    end

    
    println("Collecting results")
    pz = exp.(-0.5*chi2grid)
    cpz = cumsum(pz, dims=2) ./ sum(pz, dims=2)
    z_l95 = zgrid[map(argmin, eachrow(abs.(cpz .- 0.025)))]
    z_l68 = zgrid[map(argmin, eachrow(abs.(cpz .- 0.160)))]
    z_med = zgrid[map(argmin, eachrow(abs.(cpz .- 0.500)))]
    z_u68 = zgrid[map(argmin, eachrow(abs.(cpz .- 0.840)))]
    z_u95 = zgrid[map(argmin, eachrow(abs.(cpz .- 0.975)))]
    
    izbest = map(argmin, eachrow(chi2grid))
    zbest = zgrid[izbest]
    coeffsbest = zeros(nobj,ntempl)
    @threads for j in 1:nobj
        coeffsbest[j,:] = coeffs[j,izbest[j],:]
    end
    chi2best = vec(minimum(chi2grid, dims=2))
    
    bad_objs = sum(chi2grid .== -1, dims=2) .== nz
    zbest[bad_objs] .= -1
    chi2best[bad_objs] .= -1

    println("Writing summary output to $output_file")
    write_data(output_file, 
               ["ID", "z_best", "chi2", "z_l95", "z_l68", "z_med", "z_u68", "z_u95"], 
               Any[IDs, zbest, chi2best, z_l95, z_l68, z_med, z_u68, z_u95])

    if output_pz
        println("Writing P(z) output to $output_file_pz")
        temp_pz = pz ./ trapz(zgrid, pz)
        temp_pz = vcat(transpose(zgrid), temp_pz)
        temp_IDs = vcat([-1], IDs)
        write_data(output_file_pz, 
                   ["ID", "Pz"], 
                   Any[temp_IDs, temp_pz])
    end

    if output_templ
        println("Writing template output to $output_file_templ")

        fnu_templ = zeros(nobj, length(templwav))
        idx  = searchsortedfirst(templwav, 1215.67)
        idx_igm = searchsortedfirst(igm_wavelengths, 1215.67)

        @threads for j in 1:nobj
            z = zbest[j]
            i = argmin(abs.(zgrid .- z))

            wav_obs = templwav .* (1+z)
            
            # Interpolate the IGM transmission at this redshift
            iz_up = searchsortedfirst(igm_redshifts, z)
            transmission = igm_transmission[iz_up,:]
            interp = linear_interpolation([0.0;igm_wavelengths[1:idx_igm-1]], [0.0;transmission[1:idx_igm-1]], extrapolation_bc=Flat())
            y1 = interp(templwav[1:idx-1])
            interp = linear_interpolation([igm_wavelengths[idx_igm:end];1226.0], [transmission[idx_igm:end];1.0], extrapolation_bc=Flat())
            y2 = interp(templwav[idx:end])
            transmission  = [y1; y2]
            
            fnu_j = templfnu * coeffsbest[j,:]
            fnu_j = fnu_j .* transmission
            fnu_templ[j,:] = fnu_j

        end
        temp_templ = vcat(transpose(templwav), fnu_templ)
        temp_IDs = vcat([-1], IDs)
        temp_zbest = vcat([-1], zbest)
        write_data(output_file_templ, 
                   ["ID", "z_best", "templ"], 
                   Any[temp_IDs, temp_zbest, temp_templ])
    end
                
    


    return errno
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
    # filter_files = glob("*.dat", filterpath)
    println("Available filters:")

    for filter_name in filter_names
        filter_description = filter_directory[filter_name]["description"]
        print(rpad(filter_name, 25))
        println(filter_description)
        if !isfile(filterpath * filter_name)
            panic("Filter file $filter_name not found.")
        end
    end
end

function get_filter(nickname)
    """
    Get the filter name from the nickname.
    """
    if haskey(filter_nicknames, nickname)
        real_name = filter_nicknames[nickname]
        filt = readdlm(filterpath * real_name)
        return filt[:,1], filt[:,2]
    else
        return panic("Filter '$nickname' not found.")
    end
end

function get_pivot_wavelengths(bands)
    """
    Compute the pivot wavelengths for the given bands.
    """
    nband = length(bands)
    pivot_wavs = zeros(nband)
    for k in 1:nband
        band = bands[k]
        fwav, ftrans = get_filter(band)
        pivot_wavs[k] = sqrt(trapz(fwav, ftrans) / trapz(fwav, ftrans ./ (fwav .^ 2)))
    end
    return pivot_wavs
end


function list_templates()
    """
    List the available templates in the template_set directory.
    """
    println("Available template sets:")
    for template_set_name in template_set_names
        template_files = template_directory[template_set_name]["files"]
        for template_file in template_files
            if !isfile(templatepath * template_file)
                panic("Template file $template_file not found.")
            end
        end
        nfiles = length(template_files)
        println("* $template_set_name ($nfiles)")
    end
end





Base.@main

