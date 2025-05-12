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
    
    # Reset errno
    global errno = 0
    
    io = stdout
    print_ascii(io)
    check_version()
    
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
 
    templates = [joinpath(templatepath, file) for file in template_directory[template_set]["files"]]
    ntempl = length(templates)
    if ntempl == 0
        panic("No templates found in the template set $template_set.")
    end


    # Print out the templates for the user 
    for (i, templ) in enumerate(templates)
        templwav, templflux, templz = load_template(templ)
        templ_shortname = basename(templ)

        if templz === nothing
            s = length(templflux)
            println("Template $i: $templ_shortname ($s,)")
        else
            s = size(templflux)
            println("Template $i: $templ_shortname ($s)")
        end
    end
    
    # Build template grid 
    templgrid = zeros(ntempl, nz, nband)
    println("Building template grid: " * summary(templgrid))
    
    # Load in IGM transmission data
    if !haskey(fitting, "igm")
        error("parameter `igm` not found in the parameter file.")
    end
    
    igm_file = h5open(igmpath * fitting["igm"] * ".hdf5", "r")
    igm_redshifts = igm_file["redshifts"][:]
    igm_wavelengths = igm_file["wavelengths"][:]
    igm_transmission = igm_file["transmission"][:,:]
    close(igm_file)
    idx_igm = searchsortedfirst(igm_wavelengths, 1215.67)

    iter = ProgressBar(1:ntempl)
    for i in iter

        templ_shortname = basename(templates[i])
        templwav_i, templfnu_i, templz_i = load_template(templates[i])
        idx  = searchsortedfirst(templwav_i, 1215.67)

        @threads for j in 1:nz
            z = zgrid[j]

            wav_obs = templwav_i .* (1+z)
            
            # Interpolate the IGM transmission at this redshift
            iz_up = searchsortedfirst(igm_redshifts, z)
            transmission = igm_transmission[iz_up,:]

            interp = linear_interpolation([0.0;igm_wavelengths[1:idx_igm-1]], [0.0;transmission[1:idx_igm-1]], extrapolation_bc=Flat())
            y1 = interp(templwav_i[1:idx-1])
            interp = linear_interpolation([igm_wavelengths[idx_igm:end];1226.0], [transmission[idx_igm:end];1.0], extrapolation_bc=Flat())
            y2 = interp(templwav_i[idx:end])
            transmission  = [y1; y2]

            if templz_i === nothing
                templfnu_j = templfnu_i .* transmission
            else
                zindex = argmin(abs.(templz_i .- z))
                templfnu_j = templfnu_i[:,zindex] .* transmission
            end
            
            @threads for k in 1:nband
                
                band = bands[k]
                fwav, ftrans = get_filter(band)
                nu = 1 ./ fwav
                interp = linear_interpolation(wav_obs, templfnu_j)
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

            valid = isfinite.(fnu_j) .&  isfinite.(efnu_tot_j) .& (efnu_tot_j .> 0.0)
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
            chi2grid[j,i] = sum(((fnu_j[valid] .- fnu_mod_j[valid]) .^ 2) ./ efnu_tot_j[valid] .^ 2)
            
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

    photobest = zeros(nobj, nband)
    @threads for j in 1:nobj
        photobest[j,:] = transpose(templgrid[:,izbest[j],:]) * coeffsbest[j,:]
    end
    
    bad_objs = sum(chi2grid .== -1, dims=2) .== nz
    zbest[bad_objs] .= -1
    chi2best[bad_objs] .= -1

    println("Writing summary output to $output_file")
    cols = ["ID", "z_best", "chi2", "z_l95", "z_l68", "z_med", "z_u68", "z_u95"] 
    data = Any[IDs, zbest, chi2best, z_l95, z_l68, z_med, z_u68, z_u95]
    for (i, band) in enumerate(bands)
        push!(cols, band)
        push!(data, photobest[:,i])
    end
    write_data(output_file, cols, data)

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

        # Everything will get interpolated to a common wavelength grid 
        # Defined by the template with the largest wavelength array 
        # (i.e., highest resolution/largest range)
        common_templwav = nothing
        for (i, templ) in enumerate(templates)
            templwav_i, templfnu_i, templz_i = load_template(templ)
            if common_templwav === nothing
                common_templwav = templwav_i
            elseif length(templwav_i) > length(common_templwav)
                common_templwav = templwav_i
            end
        end
        
        idx  = searchsortedfirst(common_templwav, 1215.67)
        nwav = length(common_templwav)
                
        # Build template grid 
        templgrid = zeros(ntempl, nz, nwav)
    
        igm_file = h5open(igmpath * fitting["igm"] * ".hdf5", "r")
        igm_redshifts = igm_file["redshifts"][:]
        igm_wavelengths = igm_file["wavelengths"][:]
        igm_transmission = igm_file["transmission"][:,:]
        close(igm_file)
        idx_igm = searchsortedfirst(igm_wavelengths, 1215.67)

        igm_grid = zeros(nz, nwav)
        for i in 1:nz
            z = zgrid[i]
            
            # Interpolate the IGM transmission at this redshift
            iz_up = searchsortedfirst(igm_redshifts, z)
            t = igm_transmission[iz_up,:]

            interp = linear_interpolation([0.0;igm_wavelengths[1:idx_igm-1]], [0.0;t[1:idx_igm-1]], extrapolation_bc=Flat())
            y1 = interp(common_templwav[1:idx-1])
            interp = linear_interpolation([igm_wavelengths[idx_igm:end];1226.0], [t[idx_igm:end];1.0], extrapolation_bc=Flat())
            y2 = interp(common_templwav[idx:end])
            t  = [y1; y2]

            igm_grid[i,:] = t
        end

        templfnu = zeros(nwav, nz, ntempl)
        for i in 1:ntempl

            templwav_i, templfnu_i, templz_i = load_template(templates[i])

            if templz_i === nothing
                # Interpolate the template to the common wavelength grid
                interp = linear_interpolation(templwav_i, templfnu_i, extrapolation_bc=Flat())
                templfnu_i = interp(common_templwav)
            end

            @threads for j in 1:nz
                z = zgrid[j]

                transmission = igm_grid[j,:]
                
                if templz_i === nothing
                    templfnu_j = templfnu_i .* transmission
                else
                    zindex = argmin(abs.(templz_i .- z))
                    interp = linear_interpolation(templwav_i, templfnu_i[:,zindex], extrapolation_bc=Flat())
                    templfnu_j = interp(common_templwav) .* transmission
                end
                
                templfnu[:,j,i] = templfnu_j
                    
            end
        end
    
        fnu_templ = zeros(nobj, nwav)
        iter = ProgressBar(1:nobj)
        for i in iter
            z = zbest[i]
            j = izbest[i] 
            fnu_templ[i,:] = templfnu[:,j,:] * coeffsbest[i,:]
        end

        temp_templ = vcat(transpose(common_templwav), fnu_templ)
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

function make_bins(wavs)
    """ Given a series of wavelength points, find the edges and widths
    of corresponding wavelength bins. """
    edges = zeros(length(wavs)+1)
    widths = zeros(length(wavs))
    edges[1] = wavs[1] - (wavs[2] - wavs[1])/2
    widths[end] = (wavs[end] - wavs[end-1])
    edges[end] = wavs[end] + (wavs[end] - wavs[end-1])/2
    edges[2:end-1] = (wavs[2:end] + wavs[1:end-1])/2
    widths[1:end-1] = edges[2:end-1] - edges[1:end-2]
    return edges, widths
end

function spectres(new_wavs, old_wavs, old_fluxes; old_errs=nothing, fill_value=0.0)
    """
    Interpolate a spectrum to a new wavelength grid.
    """
    # interp = linear_interpolation(wav_old, flux_old, extrapolation_bc=Flat())
    # flux_new = interp(wav_new)
    # return flux_new

    # Make arrays of edge positions and widths for the old and new bins

    old_edges, old_widths = make_bins(old_wavs)
    new_edges, new_widths = make_bins(new_wavs)

    # Generate output arrays to be populated
    new_fluxes = zeros(length(new_wavs))

    if !isnothing(old_errs)
        if length(old_errs) != length(old_fluxes)
            panic("old_fluxes and old_errs must be the same length")
        else
            new_errs = copy(new_fluxes)
        end
    end

    start = 1
    stop = 1

    # Calculate new flux and uncertainty values, looping over new bins
    for j in 1:length(new_wavs)

        # Add filler values if new_wavs extends outside of spec_wavs
        if (new_edges[j] < old_edges[1]) || (new_edges[j+1] > old_edges[end])
            new_fluxes[j] = fill_value

            if !isnothing(old_errs)
                new_errs[j] = fill_value
            end

            if (j == 0) || (j == length(new_wavs))
                # warnings.warn(
                #     "Spectres: new_wavs contains values outside the range "
                #     "in spec_wavs, new_fluxes and new_errs will be filled "
                #     "with the value set in the 'fill' keyword argument "
                #     "(by default 0).",
                #     category=RuntimeWarning,
                # )
                continue
            end
        end

        # Find first old bin which is partially covered by the new bin
        while old_edges[start+1] <= new_edges[j]
            start += 1
        end

        # Find last old bin which is partially covered by the new bin
        while old_edges[stop+1] < new_edges[j+1]
            stop += 1
        end

        # If new bin is fully inside an old bin start and stop are equal
        if stop == start
            new_fluxes[j] = old_fluxes[start]
            if !isnothing(old_errs)
                new_errs[j] = old_errs[start]
            end
            
            # Otherwise multiply the first and last old bin widths by P_ij
        else
            start_factor = ((old_edges[start+1] - new_edges[j]) / (old_edges[start+1] - old_edges[start]))
            
            end_factor = ((new_edges[j+1] - old_edges[stop]) / (old_edges[stop+1] - old_edges[stop]))
            
            old_widths[start] *= start_factor
            old_widths[stop] *= end_factor
            
            # Populate new_fluxes spectrum and uncertainty arrays
            f_widths = old_widths[start:stop+1] .* old_fluxes[start:stop+1]
            new_fluxes[j] = sum(f_widths)/sum(old_widths[start:stop+1])
            
            if !isnothing(old_errs)
                e_wid = old_widths[start:stop+1]*old_errs[start:stop+1]
                
                new_errs[j] = sqrt(sum(e_wid .^2)) / sum(old_widths[start:stop+1])
                
                # Put back the old bin widths to their initial values
                old_widths[start] /= start_factor
                old_widths[stop] /= end_factor
            end
        end
    end
                
    # If errors were supplied return both new_fluxes and new_errs.
    if !isnothing(old_errs)
        return new_fluxes, new_errs
    else
        return new_fluxes
    end
end


function load_template(file)
    """
    Load a template from a file.
    """

    isfits = endswith(file, ".fits")

    # If the template is a FITS file, read the "wave" and "flux" columns 
    # and check for redshift-dependence 
    if isfits
        t = FITS(file)
        templwav = read(t[2], "wave")
        templflam = read(t[2], "flux")

        # Check for redshift-dependence
        if ndims(templflam) == 2
            zdependent = true
        else
            zdependent = false
        end

        # If the template is redshift-dependent
        if zdependent
            # get the redshift values from the header
            hdr = read_header(t[2])
            templnz = hdr["NZ"]
            templz = zeros(templnz)
            for tmp in 1:templnz
                tmp0 = tmp-1
                templz[tmp] = hdr["Z$tmp0"]
            end
            
            # transpose the flux arrray 
            templflam = transpose(templflam)
        end
        templfnu = templflam .* templwav .^ 2

    else
        # If the template is an ASCII file, read the first two columns 
        t = readdlm(file)
        templwav = t[:,1]
        templflam = t[:,2]
        templfnu =  templflam .* templwav .^ 2
        templz = nothing
    end
    
    return templwav, templfnu, templz
end




Base.@main

