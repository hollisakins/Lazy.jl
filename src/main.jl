supports_color(io) = get(io, :color, false)

# Custom exception type for Lazy.jl
struct LazyError <: Exception
    msg::String
    LazyError(msg::String) = new(msg)
end

Base.showerror(io::IO, e::LazyError) = print(io, "LazyError: ", e.msg)

using TOML
using FITSIO
using NonNegLeastSquares, Glob
using DelimitedFiles
using Interpolations
using Trapz
using LinearAlgebra
using Base.Threads
using HDF5
using OrderedCollections

igmpath = @__DIR__() * "/igm_data/"
templatepath = @__DIR__() * "/templates/"
filterpath = @__DIR__() * "/filter_files/"

"""
    fit_single_object(j, fnu_j, efnu_j, templgrid, template_error_grid, zgrid, nphot_min, nband, ntempl, nz; interrupted_flag=nothing)

Fit a single object across all redshifts and return results.
Returns (zbest, chi2best, coeffsbest, chi2_row) where chi2_row is the full chi2 vs redshift for P(z) calculation.
"""
function fit_single_object(j::Int, fnu_j::Vector{Float64}, efnu_j::Vector{Float64}, 
                          templgrid::Array{Float64,3}, template_error_grid::Matrix{Float64}, 
                          zgrid::Vector{Float64}, nphot_min::Int, nband::Int, ntempl::Int, nz::Int;
                          interrupted_flag::Union{Nothing,Ref{Bool}}=nothing)
    
    # Wrap entire function in try-catch to handle individual object failures gracefully
    try
        # Pre-allocate working arrays to avoid repeated allocations
        templgrid_ij = Matrix{Float64}(undef, nband, ntempl)
        fnu_mod_j = Vector{Float64}(undef, nband)
        chi2_row = Vector{Float64}(undef, nz)
        
        # Track best fit for this object
        best_chi2 = Inf
        best_z_idx = 0
        best_coeffs = zeros(ntempl)
        
        # Loop over all redshifts for this object
        for i in 1:nz
            # Check for interruption if flag provided
            if interrupted_flag !== nothing && interrupted_flag[]
                # Fill remaining chi2 values with -1 and break
                chi2_row[i:end] .= -1
                break
            end
            
            templgrid_i = templgrid[:,i,:]
            tefz = template_error_grid[i, :]
            
            efnu_tot_j = sqrt.( efnu_j .^ 2 + (tefz .* max.(fnu_j, 0.0)) .^ 2 )
            
            valid = isfinite.(fnu_j) .&  isfinite.(efnu_tot_j) .& (efnu_tot_j .> 0.0)
            if sum(valid) < 2
                chi2_row[i] = -1
                continue
            end
            
            detect = valid .& ((fnu_j ./ efnu_j) .> 2.0)
            if sum(detect) < nphot_min
                chi2_row[i] = -1
                continue
            end
            
            snr_j = fnu_j ./ efnu_tot_j
            
            # Optimized: avoid transpose by using templgrid_i directly with broadcasting
            # templgrid_i is (ntempl, nband), we need (nband, ntempl) divided by efnu_tot_j
            for k in 1:nband
                for t in 1:ntempl
                    templgrid_ij[k, t] = templgrid_i[t, k] / efnu_tot_j[k]
                end
            end
            
            result = nonneg_lsq(templgrid_ij[valid,:], snr_j[valid] ; alg=:nnls)[:]
            
            # Numerical stability checks
            if !all(isfinite, result)
                #@warn "⚠️ Non-finite coefficients detected for object $j at z=$(zgrid[i]). Skipping this redshift."
                chi2_row[i] = -1
                continue
            end
            
            if all(result .≈ 0.0)
                #@warn "⚠️ All-zero coefficients detected for object $j at z=$(zgrid[i]). Poor template fit."
                chi2_row[i] = 1e10  # High chi2 to mark as poor fit
                continue
            end
            
            # Optimized: avoid transpose by manual matrix multiplication
            fill!(fnu_mod_j, 0.0)
            for k in 1:nband
                for t in 1:ntempl
                    fnu_mod_j[k] += templgrid_i[t, k] * result[t]
                end
            end
            chi2_j = sum(((fnu_j[valid] .- fnu_mod_j[valid]) .^ 2) ./ efnu_tot_j[valid] .^ 2)
            
            # Additional numerical stability check for chi2
            if !isfinite(chi2_j) || chi2_j < 0
                #@warn "⚠️ Invalid chi2 ($chi2_j) for object $j at z=$(zgrid[i]). Skipping."
                chi2_row[i] = -1
                continue
            end
            
            # Store chi2 for P(z) calculation
            chi2_row[i] = chi2_j
            
            # Update best fit if this is better
            if chi2_j < best_chi2
                best_chi2 = chi2_j
                best_z_idx = i
                best_coeffs = result
            end
        end
    
        # Determine final results
        if best_z_idx > 0
            zbest = zgrid[best_z_idx]
            chi2best = best_chi2
            coeffsbest = best_coeffs
        else
            zbest = -1.0
            chi2best = -1.0
            coeffsbest = zeros(ntempl)
        end
        
        return (zbest, chi2best, coeffsbest, chi2_row)
        
    catch e
        # Handle any unexpected errors for this individual object
        @warn "⚠️ Error fitting object $j: $e"
        
        # Return safe default values indicating failure
        zbest_failed = -1.0
        chi2best_failed = -1.0
        coeffsbest_failed = zeros(ntempl)
        chi2_row_failed = fill(-1.0, nz)
        
        return (zbest_failed, chi2best_failed, coeffsbest_failed, chi2_row_failed)
    end
end

# Custom Progress Bar Utility
mutable struct ProgressBar
    total::Int
    current::Atomic{Int}
    description::String
    start_time::Float64
    last_update_time::Ref{Float64}
    show_rate::Bool
    show_eta::Bool
    bar_width::Int
    prefix_emoji::String
    
    function ProgressBar(total::Int, description::String; 
                        show_rate::Bool=true, show_eta::Bool=true, 
                        bar_width::Int=20, prefix_emoji::String="🧠")
        new(total, Atomic{Int}(0), description, time(), Ref(time()), 
            show_rate, show_eta, bar_width, prefix_emoji)
    end
end

"""
Format numbers with k/M suffixes for better readability
"""
function format_number(n::Real)::String
    if n >= 1_000_000
        return "$(round(n/1_000_000, digits=1))M"
    elseif n >= 1_000
        return "$(round(n/1_000, digits=1))k"
    else
        return string(round(Int, n))
    end
end

"""
Format elapsed time in a human-readable format
"""
function format_time(elapsed::Float64)::String
    if elapsed < 1
        return "$(round(Int, elapsed * 1000))ms"
    elseif elapsed < 60
        return "$(round(elapsed, digits=1))s"
    elseif elapsed < 3600
        minutes = floor(Int, elapsed / 60)
        seconds = round(Int, elapsed % 60)
        return "$(minutes)m$(seconds)s"
    else
        hours = floor(Int, elapsed / 3600)
        minutes = floor(Int, (elapsed % 3600) / 60)
        return "$(hours)h$(minutes)m"
    end
end

"""
Update progress bar display
"""
function update!(pb::ProgressBar, new_value::Int; force::Bool=false)
    pb.current[] = new_value
    current_time = time()
    
    # Only update display if forced or enough time has passed (throttling)
    if force || (current_time - pb.last_update_time[]) > 0.1
        pb.last_update_time[] = current_time
        display_progress(pb)
    end
end

"""
Increment progress by 1
"""
function increment!(pb::ProgressBar; force::Bool=false)
    new_value = atomic_add!(pb.current, 1)
    current_time = time()
    
    # Only update display if forced or enough time has passed
    if force || (current_time - pb.last_update_time[]) > 0.1
        pb.last_update_time[] = current_time
        display_progress(pb)
    end
    
    return new_value
end

"""
Display the progress bar
"""
function display_progress(pb::ProgressBar)
    current = pb.current[]
    percentage = min(100, round(Int, (current / pb.total) * 100))
    
    # Create progress bar
    filled = round(Int, (current / pb.total) * pb.bar_width)
    bar = "█"^filled * "░"^(pb.bar_width - filled)
    
    # Format counts
    current_str = format_number(current)
    total_str = format_number(pb.total)
    
    # Build display string
    display_str = "$(pb.prefix_emoji) $(pb.description) [$bar] $percentage% ($current_str/$total_str"
    
    # Add rate and ETA if requested
    if pb.show_rate || pb.show_eta
        elapsed = time() - pb.start_time
        if elapsed > 0 && current > 0
            rate = current / elapsed
            
            if pb.show_rate
                display_str *= "; $(round(rate, digits=1)) obj/s"
            end
            
            if pb.show_eta && current < pb.total
                remaining = pb.total - current
                eta = remaining / rate
                
                eta_str = if !isfinite(eta) || eta <= 0
                    "∞"
                elseif eta < 60
                    "$(round(Int, eta))s"
                elseif eta < 3600
                    "$(round(Int, eta ÷ 60))m $(round(Int, eta % 60))s"
                else
                    "$(round(Int, eta ÷ 3600))h $(round(Int, (eta % 3600) ÷ 60))m"
                end
                
                display_str *= ", ETA: $eta_str"
            end
        end
    end
    
    display_str *= ")"
    
    # Clear line and print
    print("\r$display_str")
    flush(stdout)
end

"""
Finish the progress bar with completion message
"""
function finish!(pb::ProgressBar; final_message::String="")
    pb.current[] = pb.total
    current_time = time()
    elapsed = current_time - pb.start_time
    
    # Use elapsed time as default final message
    if final_message == ""
        final_message = "took $(format_time(elapsed))"
    end
    
    # Build final display string respecting original settings
    current_str = format_number(pb.total)
    total_str = format_number(pb.total)
    
    display_str = "$(pb.prefix_emoji) $(pb.description) [$("█"^pb.bar_width)] 100% ($current_str/$total_str"
    
    # Add rate only if show_rate was enabled
    if pb.show_rate && elapsed > 0
        rate = pb.total / elapsed
        display_str *= "; $(round(rate, digits=1)) obj/s"
    end
    
    display_str *= ", $final_message)"
    
    print("\r$display_str")
    println()  # New line after completion
    flush(stdout)
end

# Spinner utilities for dynamic feedback
const SPINNER_CHARS = ['⠋', '⠙', '⠹', '⠸', '⠼', '⠴', '⠦', '⠧', '⠇', '⠏']

function with_spinner(operation::Function, message::String, success_emoji::String="✅")
    """
    Run an operation with a spinner for visual feedback.
    Falls back to static messages in non-interactive terminals.
    """
    # Check for terminal capabilities - TTY detection doesn't work reliably through module invocation
    term_val = get(ENV, "TERM", "")
    has_proper_term = term_val != "dumb" && term_val != ""
    # Allow spinners if we have a proper TERM, even if TTY detection fails
    if has_proper_term
        # Use spinner for interactive terminals
        spinner_active = Ref(true)
        max_line_length = Ref(0)
        
        # Start spinner task
        spinner_task = @async begin
            i = 1
            while spinner_active[]
                spinner_line = "$(SPINNER_CHARS[i]) $message..."
                print("\r$spinner_line")
                max_line_length[] = max(max_line_length[], length(spinner_line))
                flush(stdout)
                sleep(0.1)
                i = (i % length(SPINNER_CHARS)) + 1
            end
        end
        
        try
            result = operation()
            spinner_active[] = false
            wait(spinner_task)
            # Clear the spinner line completely and print success message
            clear_line = "\r" * " "^max_line_length[] * "\r"
            print("$clear_line$success_emoji $message\n")
            flush(stdout)
            return result
        catch e
            spinner_active[] = false
            wait(spinner_task)
            # Clear the spinner line completely and print error message
            clear_line = "\r" * " "^max_line_length[] * "\r"
            print("$clear_line❌ $message (failed)\n")
            flush(stdout)
            rethrow(e)
        end
    else
        # Fallback for non-interactive terminals - only print completion message
        result = operation()
        println("$success_emoji $message")
        return result
    end
end

function print_help(io, cmd::String = "main")
    printstyled(io, "Lazy.jl \n", bold = true)
    if cmd == "main"
        println("usage: lazy <command> <options>")
        println("")
        println("  fit <options>       Fit the data")
        println("  list-templates      List the available templates")
        println("  list-filters        List the available filters")
        println("  cache-clear         Clear cached template grids")
        println("  -v, --version       Show version information")
        println("  -h, --help          Show this help message")
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

function load_templates()
    """
    Load the templates from the template directory.
    """
    if !isdir(templatepath)
        throw(LazyError("Template directory $templatepath not found"))
    end
    return TOML.parsefile(templatepath * "template_directory.toml")
end

function load_filters()
    """
    Load the filters from the filter directory.
    """
    if !isdir(filterpath)
        throw(LazyError("Filter directory $filterpath not found"))
    end
    return TOML.parsefile(filterpath * "filter_directory.toml")
end

function main(argv::Vector{String})
    
    io = stdout    
    if argv == []
        print_ascii(io)
        check_version()
        print_help(io)
        return 0
    end

    # Parse command line arguments
    while length(argv) > 0
        x = popfirst!(argv)
        
        # Return the help message and exit
        if x == "-h" || x == "--help"
            print_ascii(io)
            check_version()        
            print_help(io)
            return 0
        
        # Return version information
        elseif x == "-v" || x == "--version"
            println("Lazy.jl version $version")
            return 0
        
        # 
        elseif x == "fit"
            print_ascii(io)
            check_version()

            if argv == []
                print_help(io, "fit")
                return 0
            end
        
            y = popfirst!(argv)

            if y == "-p" || y == "--param"
                if length(argv) < 1
                    throw(LazyError("expected parameter file argument after `-p`"))
                end
                param = popfirst!(argv)
                return fit(param)
            elseif y == "-h" || y == "--help"
                print_help(io, "fit")
                return 0
            else
                throw(LazyError("unrecognized argument: $y"))
            end


        elseif x == "list-templates"
            # List the available templates
            print_ascii(io)
            check_version()
            return list_templates()
            
        elseif x == "list-filters"
            # List the available filters
            print_ascii(io)
            check_version()
            return list_filters()
            
        elseif x == "params"
            return params()
            
        elseif x == "cache-clear"
            # Clear template cache
            print_ascii(io)
            check_version()
            removed_count, freed_gb = clear_cache()
            if removed_count > 0
                println("Cleared $removed_count cached template grids ($(round(freed_gb, digits=2)) GB freed)")
            else
                println("No cached template grids found")
            end
            return 0
            
        else
            # Argument not recognized
            throw(LazyError("unrecognized argument: $x"))
        end
    end
    if param == nothing
        throw(LazyError("parameter file not specified. Use -p <param_file>"))
    end

end

function precompute_template_errors(zgrid::Vector{Float64}, pivot_wavs::Vector{Float64}, 
                                   tef_x::Vector{Float64}, tef_y::Vector{Float64},
                                   tef_xmin::Float64, tef_xmax::Float64, 
                                   tef_clip_min::Float64, tef_clip_max::Float64,
                                   template_error_scale::Float64)::Matrix{Float64}
    """
    Pre-compute template error grid for all redshifts and bands.
    This avoids repeated computation during by-object fitting.
    
    Returns: template_error_grid[nz, nband]
    """
    nz = length(zgrid)
    nband = length(pivot_wavs)
    template_error_grid = zeros(nz, nband)
    
    for i in 1:nz
        z = zgrid[i]
        pivot_wavs_rest = pivot_wavs ./ (1 + z)
        
        # Interpolate template error function
        tef_interp = linear_interpolation(tef_x, tef_y, extrapolation_bc=Flat())
        tefz = tef_interp(pivot_wavs_rest)
        
        # Apply clipping
        tefz[pivot_wavs_rest .< tef_xmin] .= tef_clip_min
        tefz[pivot_wavs_rest .> tef_xmax] .= tef_clip_max
        tefz = tefz * template_error_scale
        
        template_error_grid[i, :] = tefz
    end
    
    return template_error_grid
end

function estimate_memory_usage(nobj::Int, nz::Int, nband::Int, ntempl::Int)::Dict{String, Float64}
    """
    Estimate memory usage for the fitting process in GB.
    
    Returns a dictionary with estimated memory usage for different components.
    """
    # Template grid: ntempl × nz × nband × 8 bytes (Float64)
    templgrid_gb = (ntempl * nz * nband * 8) / (1024^3)
    
    # Chi2 grid: nobj × nz × 8 bytes (Float64)  
    chi2grid_gb = (nobj * nz * 8) / (1024^3)
    
    # Template error grid: nz × nband × 8 bytes (Float64)
    template_error_gb = (nz * nband * 8) / (1024^3)
    
    # Best-fit coefficients: nobj × ntempl × 8 bytes (Float64)
    coeffs_gb = (nobj * ntempl * 8) / (1024^3)
    
    # Photometric data: nobj × nband × 2 (flux + error) × 8 bytes
    photdata_gb = (nobj * nband * 2 * 8) / (1024^3)
    
    # Working arrays per thread (templgrid_ij, fnu_mod_j per thread)
    nthreads = Threads.nthreads()
    working_arrays_gb = (nthreads * (nband * ntempl + nband) * 8) / (1024^3)
    
    # Total estimated peak usage
    peak_gb = templgrid_gb + chi2grid_gb + template_error_gb + coeffs_gb + photdata_gb + working_arrays_gb
    
    return Dict(
        "template_grid" => templgrid_gb,
        "chi2_grid" => chi2grid_gb, 
        "template_error_grid" => template_error_gb,
        "coefficients" => coeffs_gb,
        "photometric_data" => photdata_gb,
        "working_arrays" => working_arrays_gb,
        "estimated_peak" => peak_gb
    )
end

function print_memory_estimate(mem_dict::Dict{String, Float64})
    """
    Print a formatted memory usage estimate.
    """
    println("====================================")
    println("Memory Usage Estimate:")
    println("  Template grid:      $(round(mem_dict["template_grid"], digits=2)) GB")
    println("  Chi2 grid:          $(round(mem_dict["chi2_grid"], digits=2)) GB") 
    println("  Template errors:    $(round(mem_dict["template_error_grid"], digits=2)) GB")
    println("  Coefficients:       $(round(mem_dict["coefficients"], digits=2)) GB")
    println("  Photometric data:   $(round(mem_dict["photometric_data"], digits=2)) GB")
    println("  Working arrays:     $(round(mem_dict["working_arrays"], digits=2)) GB")
    println("  ----------------------------------------")
    println("  Estimated peak:     $(round(mem_dict["estimated_peak"], digits=2)) GB")
    
    if mem_dict["estimated_peak"] > 16.0
        println("  ⚠️  WARNING: High memory usage detected!")
        println("     Consider reducing parameter space if system has <$(round(mem_dict["estimated_peak"] * 1.2, digits=1))GB RAM")
    elseif mem_dict["estimated_peak"] > 8.0
        println("  ⚠️  CAUTION: Moderate memory usage.")
        println("     Monitor system memory during fitting.")
    else
        println("  ✅ Memory usage looks reasonable.")
    end
    println("====================================")
end

function fit(param)

    # Load in TOML parameter file
    if !isfile(param)
        throw(LazyError("parameter file $param not found"))
    end
    if !endswith(param, ".toml")
        throw(LazyError("parameter file $param must have a .toml extension"))
    end
    println("Loading parameter file: $param")
    param = TOML.parsefile(param)
    

    # Attempt to read in the input catalog file
    if haskey(param, "io")
        io = param["io"]
        if haskey(io, "input_catalog")
            input_catalog = io["input_catalog"]
            cat = with_spinner(() -> FITS(input_catalog), "Loading catalog: $input_catalog", "📊")
        else
            throw(LazyError("parameter `io.input_catalog` not found in the parameter file"))
        end

        if haskey(io, "output_file")
            output_file = io["output_file"]
        else
            throw(LazyError("parameter `io.output_file` not found in the parameter file"))
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
    else
        throw(LazyError("section `io` not found in the parameter file"))
    end

    IDs = Int.(read(cat[2], "ID"))
    nobj = length(IDs)
    println("📊 Objects: $nobj")
        
    println("====================================")
    if !haskey(param, "translate")
        throw(LazyError("section `translate` not found in the parameter file"))
    end
    
    translate = param["translate"]
    bands = sort(collect(keys(param["translate"])))
    # Load in the data
    fnu, efnu, bands = with_spinner(() -> load_data(cat, bands, translate), "Loading photometric data", "📊")
    nband = length(bands)
    println("📊 Bands: $nband")
    println("📊 Filters: $bands")

    println("📊 Flux data: " * summary(fnu))
    println("📊 Error data: " * summary(efnu))

    if !haskey(io, "missing_data_format")
        throw(LazyError("parameter `io.missing_data_format` not found in the parameter file"))
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
        println("⚠️ Unrecognized missing data format: $missing_data_format, using input catalog as is")
    end

        


    println("====================================")

    
    if !haskey(param, "fitting")
        throw(LazyError("section `fitting` not found in the parameter file"))
    end

    fitting = param["fitting"]

    if !haskey(fitting, "nphot_min")
        throw(LazyError("parameter `fitting.nphot_min` not found in the parameter file"))
    end
    nphot_min = fitting["nphot_min"]

    if !haskey(fitting, "sys_err")
        throw(LazyError("parameter `fitting.sys_err` not found in the parameter file"))
    end
    sys_err = fitting["sys_err"]
    fnu, efnu = set_sys_err(fnu, efnu, sys_err)

    if !haskey(fitting, "z_min")
        throw(LazyError("parameter `fitting.z_min` not found in the parameter file"))
    end
    if !haskey(fitting, "z_max")
        throw(LazyError("parameter `fitting.z_max` not found in the parameter file"))
    end
    if !haskey(fitting, "z_step")
        throw(LazyError("parameter `fitting.z_step` not found in the parameter file"))
    end

    z_min = fitting["z_min"]
    z_max = fitting["z_max"]
    z_step = fitting["z_step"]
    println("🎆 Redshift range: $z_min - $z_max (step: $z_step)")

    zgrid = collect(range(z_min, stop=z_max, step=z_step))
    nz = length(zgrid)

    println("====================================")

    if !haskey(fitting, "template_set")
        throw(LazyError("parameter `template_set` not found in the parameter file"))
    end

    template_set = fitting["template_set"]
    
    if template_set isa String
        template_directory = with_spinner(() -> load_templates(), "Loading template directory", "⚙️")
        if !haskey(template_directory, template_set)
            throw(LazyError("Template set $template_set not found in the template directory"))
        end
        println("⚙️ Template set: $template_set")
        templates = [joinpath(templatepath, file) for file in template_directory[template_set]["files"]]
    elseif template_set isa Vector{String}
        templates = [joinpath(templatepath, file) for file in template_set]
    else
        throw(LazyError("parameter `template_set` must be a string, defining the template set name, or a vector of strings, defining the paths to the template files to use"))
    end

    ntempl = length(templates)
    for templ in templates
        if !isfile(templ)
            throw(LazyError("Template file $templ not found in the templates directory"))
        end
    end

    # Estimate and report memory usage
    memory_estimate = estimate_memory_usage(nobj, nz, nband, ntempl)
    print_memory_estimate(memory_estimate)

    # Print out the templates for the user 
    for (i, templ) in enumerate(templates)
        templwav, templflux, templz = load_template(templ)
        templ_shortname = basename(templ)

        if templz === nothing
            s = length(templflux)
            println("⚙️ Template $i: $templ_shortname ($s,)")
        else
            s = size(templflux)
            println("⚙️ Template $i: $templ_shortname ($s)")
        end
    end
    
    # Check for template cache before building grid
    use_cache = get(fitting, "template_cache", true)  # Default enabled
    templgrid = nothing
    template_error_grid = nothing
    cache_key = ""
    
    if use_cache
        # Generate cache key from all relevant parameters
        cache_key = generate_cache_key(templates, zgrid, bands, fitting["igm_model"], 
                                     fitting["template_error"], fitting["template_error_scale"])
        
        # Try to load from cache
        cached_data = load_template_cache(cache_key)
        if cached_data !== nothing
            templgrid, template_error_grid = cached_data
            println("✅ Loading cached template grid ($(cache_key)) - " * summary(templgrid))
        end
    end
    
    # Build template grid if not loaded from cache
    if templgrid === nothing
        templgrid = zeros(ntempl, nz, nband)
        templgrid_size_gb = sizeof(templgrid) / (1024^3)
        println("🧠 Building template grid: " * summary(templgrid) * " (~$(round(templgrid_size_gb, digits=2)) GB)")
        if templgrid_size_gb > 8.0
            println("⚠️ Large template grid detected. Consider reducing nz, nband, or ntempl if memory issues occur.")
        end
        
        # Initialize progress bar for template processing
        template_progress = ProgressBar(ntempl, "Loading templates...", 
                                       show_rate=false, show_eta=false, 
                                       prefix_emoji="📊")
        
        # Load in IGM transmission data
    if !haskey(fitting, "igm_model")
        throw(LazyError("parameter `igm_model` not found in the parameter file"))
    end
    
    igm_redshifts, igm_wavelengths, igm_transmission = with_spinner(() -> begin
        igm_file = h5open(igmpath * fitting["igm_model"] * ".hdf5", "r")
        redshifts = igm_file["redshifts"][:]
        wavelengths = igm_file["wavelengths"][:]
        transmission = igm_file["transmission"][:,:]
        close(igm_file)
        return redshifts, wavelengths, transmission
    end, "Loading IGM model: $(fitting["igm_model"])", "🌌")
    idx_igm = searchsortedfirst(igm_wavelengths, 1215.67)
    maxz = igm_redshifts[end]
    pr = true

    for i in 1:ntempl

        templ_shortname = basename(templates[i])
        templwav_i, templfnu_i, templz_i = load_template(templates[i])
        idx  = searchsortedfirst(templwav_i, 1215.67)
        
        @threads for j in 1:nz
            z = zgrid[j]

            wav_obs = templwav_i .* (1+z)
            
            # Interpolate the IGM transmission at this redshift
            if z > maxz
                if pr
                    println("⚠️ Redshift $z is greater than the maximum redshift in the IGM model. Using IGM transmission at z=$maxz")
                    pr = false
                end
                iz_up = length(igm_redshifts)
            else
                iz_up = searchsortedfirst(igm_redshifts, z)
            end

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

            # normalize template to rest-1micron (helps with numerical stability)
            windex = argmin(abs.(templwav_i .- 1e4))
            templfnu_j /= templfnu_j[windex]

            # Sequential loop over bands (removed nested threading for thread safety)
            for k in 1:nband
                
                band = bands[k]
                fwav, ftrans = get_filter(band)
                nu = 1 ./ fwav
                interp = linear_interpolation(wav_obs, templfnu_j)
                fnu_interp = interp(fwav)
            
                result = trapz(nu, fnu_interp .* ftrans ./ nu) / trapz(nu, ftrans ./ nu)
                templgrid[i,j,k] = result

            end
        end
        
        # Update template progress
        increment!(template_progress, force=true)
    end
        # Finish template progress bar
        finish!(template_progress)
    end  # End of if templgrid === nothing
    
    # Build template error grid if not loaded from cache
    if template_error_grid === nothing
        if !haskey(fitting, "template_error")
            throw(LazyError("parameter `template_error` not found in the parameter file"))
        end
       
        if !haskey(fitting, "template_error_scale")
            throw(LazyError("parameter `template_error_scale` not found in the parameter file"))
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

        println("Building template error grid..")
        template_error_grid = precompute_template_errors(zgrid, pivot_wavs, tef_x, tef_y, 
            tef_xmin, tef_xmax, tef_clip_min, tef_clip_max, template_error_scale)
        
        # Save to cache if grid was just built and caching is enabled
        if use_cache && cache_key != ""
            metadata = Dict(
                "templates" => templates,
                "zgrid" => zgrid,
                "bands" => bands,
                "igm_model" => fitting["igm_model"],
                "template_error" => fitting["template_error"],
                "template_error_scale" => fitting["template_error_scale"],
                "ntempl" => ntempl,
                "nz" => nz,
                "nband" => nband
            )
            
            save_template_cache(templgrid, template_error_grid, metadata, cache_key)
            println("💾 Template grid cached for future use ($(cache_key))")
        end
    end

    println("===========================")
    
    # Template grids are now ready (either from cache or newly built)
    
    println("🧠 Fitting by object")
    # Initialize arrays - keep chi2grid for P(z), eliminate massive coeffs array
    chi2grid = zeros(nobj, nz)  # Keep for P(z) - only (nobj × nz)
    zbest = zeros(nobj)
    chi2best = zeros(nobj) 
    coeffsbest = zeros(nobj, ntempl)  # Only store best-fit coefficients
    
    # Memory usage warning
    chi2grid_size_gb = sizeof(chi2grid) / (1024^3)
    total_fitting_memory_gb = chi2grid_size_gb + sizeof(coeffsbest) / (1024^3)
    if total_fitting_memory_gb > 4.0
        println("⚠️ Large fitting arrays detected (~$(round(total_fitting_memory_gb, digits=2)) GB). Monitor memory usage.")
    end
    
    # Process each object individually using @spawn for better control
    # println("🧠 Fitting $nobj objects (press Ctrl+C to interrupt)")
    
    # Set up interrupt handling
    interrupted = Ref(false)
    
    # Set up signal handler for clean interruption
    try
        # Set up interrupt detection (simplified approach)
        
        # Launch fitting tasks with dynamic work distribution
        nthreads_to_use = min(Threads.nthreads(), nobj)
        
        # Optimize batch size: balance between task overhead and load balancing
        # Use larger batches for many objects, smaller for few objects
        objects_per_task = max(1, min(100, nobj ÷ (nthreads_to_use * 4)))
        ntasks = cld(nobj, objects_per_task)  # ceiling division
        
        println("⚙️ Using $nthreads_to_use threads, $ntasks tasks ($objects_per_task objects/task)")
        
        # Initialize progress bar for fitting
        progress_bar = ProgressBar(nobj, "Fitting...", 
                                  show_rate=true, show_eta=true, 
                                  prefix_emoji="🧠")
        
        # Create batched tasks for better efficiency
        tasks = Vector{Task}(undef, ntasks)
        for task_id in 1:ntasks
            start_obj = (task_id - 1) * objects_per_task + 1
            end_obj = min(task_id * objects_per_task, nobj)
            object_range = start_obj:end_obj
            
            tasks[task_id] = Threads.@spawn begin
                # Results for this batch
                batch_results = Vector{Tuple{Int, Float64, Float64, Vector{Float64}, Vector{Float64}}}()
                
                for j in object_range
                    # Check for interruption before starting work
                    if interrupted[]
                        push!(batch_results, (j, -1.0, -1.0, zeros(ntempl), fill(-1.0, nz)))
                        continue
                    end
                    
                    # Fit this object
                    zbest_j, chi2best_j, coeffsbest_j, chi2_row_j = fit_single_object(
                        j, fnu[j,:], efnu[j,:], templgrid, template_error_grid, 
                        zgrid, nphot_min, nband, ntempl, nz; interrupted_flag=interrupted
                    )
                    
                    push!(batch_results, (j, zbest_j, chi2best_j, coeffsbest_j, chi2_row_j))
                    
                    # Update progress bar
                    increment!(progress_bar)
                end
                
                return batch_results
            end
        end
        
        # Collect results as they complete
        for (task_id, task) in enumerate(tasks)
            if interrupted[] && !istaskdone(task)
                # Skip waiting for unfinished tasks if interrupted
                continue
            end
            
            # Since individual object errors are now handled within fit_single_object,
            # we can safely fetch results without batch-level error handling
            batch_results = fetch(task)
            
            # Store results from this batch
            for (j, zbest_j, chi2best_j, coeffsbest_j, chi2_row_j) in batch_results
                zbest[j] = zbest_j
                chi2best[j] = chi2best_j
                coeffsbest[j,:] = coeffsbest_j
                chi2grid[j,:] = chi2_row_j
            end
        end
        
        # Finish progress bar
        if interrupted[]
            finish!(progress_bar, final_message="interrupted")
        else
            finish!(progress_bar)
        end
        
    catch e
        if isa(e, InterruptException)
            interrupted[] = true
            println("\n⚠️ Fitting interrupted by user. Processing results for completed objects...")
        else
            rethrow(e)
        end
    end

    
    pz, z_l95, z_l68, z_med, z_u68, z_u95 = with_spinner(() -> begin
        pz = exp.(-0.5*chi2grid)
        cpz = cumsum(pz, dims=2) ./ sum(pz, dims=2)
        z_l95 = zgrid[map(argmin, eachrow(abs.(cpz .- 0.025)))]
        z_l68 = zgrid[map(argmin, eachrow(abs.(cpz .- 0.160)))]
        z_med = zgrid[map(argmin, eachrow(abs.(cpz .- 0.500)))]
        z_u68 = zgrid[map(argmin, eachrow(abs.(cpz .- 0.840)))]
        z_u95 = zgrid[map(argmin, eachrow(abs.(cpz .- 0.975)))]
        return pz, z_l95, z_l68, z_med, z_u68, z_u95
    end, "Calculating redshift statistics", "✅")
    
    bad_objs = sum(chi2grid .== -1, dims=2) .== nz
    zbest[bad_objs] .= -1
    chi2best[bad_objs] .= -1

    photobest = with_spinner(() -> begin
        photobest = zeros(nobj, nband)
        @threads for j in 1:nobj
            if zbest[j] > 0  # Only calculate for valid fits
                z_idx = argmin(abs.(zgrid .- zbest[j]))
                # Optimized: avoid transpose by manual matrix multiplication
                for k in 1:nband
                    photobest[j, k] = 0.0
                    for t in 1:ntempl
                        photobest[j, k] += templgrid[t, z_idx, k] * coeffsbest[j, t]
                    end
                end
            end
        end
        return photobest
    end, "Calculating best-fit model photometry", "✅")

    with_spinner(() -> begin
        data = OrderedDict{String, Dict{String, Any}}()
        data["ID"] = Dict("format" => "K", "data" => IDs)
        data["z_best"] = Dict("format" => "E", "data" => zbest)
        data["chi2"] = Dict("format" => "E", "data" => chi2best)
        data["z_l95"] = Dict("format" => "E", "data" => z_l95)
        data["z_l68"] = Dict("format" => "E", "data" => z_l68)
        data["z_med"] = Dict("format" => "E", "data" => z_med)
        data["z_u68"] = Dict("format" => "E", "data" => z_u68)
        data["z_u95"] = Dict("format" => "E", "data" => z_u95)
        for (i, band) in enumerate(bands)
            data[band] = Dict("format" => "E", "unit" => "fnu", "data" => photobest[:,i])
        end
        data["coeffs"] = Dict("format" => "$(ntempl)E", "data" => coeffsbest)

        write_data(output_file, data, "SUMMARY")
    end, "Writing summary output to $output_file:SUMMARY", "💾")

    if output_pz
        with_spinner(() -> begin
            temp_pz = pz ./ trapz(zgrid, pz)
            temp_pz = vcat(transpose(zgrid), temp_pz)
            temp_IDs = vcat([-1], IDs)
            data = OrderedDict{String, Dict{String, Any}}()
            data["ID"] = Dict("format" => "K", "data" => temp_IDs)
            data["Pz"] = Dict("format" => "$(nz)E", "data" => temp_pz)
            write_data(output_file, data, "PZ")
        end, "Writing P(z) output to $output_file:PZ", "💾")
    end

    if output_templ
        GC.gc()
        with_spinner(() -> begin
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
                    
        
            igm_file = h5open(igmpath * fitting["igm_model"] * ".hdf5", "r")
            igm_redshifts = igm_file["redshifts"][:]
            igm_wavelengths = igm_file["wavelengths"][:]
            igm_transmission = igm_file["transmission"][:,:]
            close(igm_file)
            idx_igm = searchsortedfirst(igm_wavelengths, 1215.67)

            igm_grid = zeros(nz, nwav)
            maxz = igm_redshifts[end]
            pr = true

            for i in 1:nz
                z = zgrid[i]
                

                # Interpolate the IGM transmission at this redshift
                if z > maxz
                    if pr
                        println("⚠️ Redshift $z is greater than the maximum redshift in the IGM model. Using IGM transmission at z=$maxz")
                        pr = false
                    end
                    iz_up = length(igm_redshifts)
                else
                    iz_up = searchsortedfirst(igm_redshifts, z)
                end
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

                for j in 1:nz
                    z = zgrid[j]

                    transmission = igm_grid[j,:]
                    
                    if templz_i === nothing
                        templfnu_j = templfnu_i .* transmission
                    else
                        zindex = argmin(abs.(templz_i .- z))
                        interp = linear_interpolation(templwav_i, templfnu_i[:,zindex], extrapolation_bc=Flat())
                        templfnu_j = interp(common_templwav) .* transmission
                    end

                    windex = argmin(abs.(common_templwav .- 1e4))
                    templfnu_j /= templfnu_j[windex]
                    
                    templfnu[:,j,i] = templfnu_j
                        
                end
            end
            
            data = OrderedDict{String, Dict{String, Any}}()
            temp_zgrid = vcat([-1], zgrid)
            data["z"] = Dict("format" => "E", "data" => temp_zgrid)
            for (i, templ) in enumerate(templates)
                templ_shortname = basename(templ)
                templ_shortname = replace(templ_shortname, ".fits" => "")
                temp_templfnu = vcat(transpose(common_templwav), transpose(templfnu[:,:,i]))
                data[templ_shortname] = Dict("format" => "$(nwav)E", "data" => temp_templfnu)
            end
            write_data(output_file, data, "TEMPL")
        end, "Writing template output to $output_file:TEMPL", "💾")
    
    end
                
    


    return 0
end

    

function load_data(cat, bands::Vector{String}, translate::Dict)::Tuple{Matrix{Float64}, Matrix{Float64}, Vector{String}}
    """
    Load the photometric data from the catalog.
    
    Returns (fnu, efnu, bands) where fnu and efnu are nobj×nband matrices.
    """
    IDs = read(cat[2], "ID")
    nobj = length(IDs)

    # Downselect to bands that are present in the catalog
    all_cols = FITSIO.colnames(cat[2])
    bands = filter(b -> (translate[b]["flux"] in all_cols) && (translate[b]["error"] in all_cols), bands)
    nband = length(bands)

    if nband == 0
        throw(LazyError("No valid bands found in the catalog. Check the translate section in the parameter file."))
    end

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
    return fnu, efnu, bands
end


function set_sys_err(fnu::Matrix{Float64}, efnu::Matrix{Float64}, sys_err::Float64)::Tuple{Matrix{Float64}, Matrix{Float64}}
    """
    Set the systematic error on the fluxes.
    
    Returns (fnu, efnu) with systematic errors added in quadrature.
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
    filter_directory = load_filters()
    filter_names = sort(collect(keys(filter_directory)))
    
    println("Available filters:")
    for filter_name in filter_names
        filter_description = filter_directory[filter_name]["description"]
        print(rpad(filter_name, 25))
        println(filter_description)
        if !isfile(filterpath * filter_name)
            throw(LazyError("Filter file $filter_name not found"))
        end
    end
end

function get_filter(nickname::String)::Tuple{Vector{Float64}, Vector{Float64}}
    """
    Get the filter transmission curve from the nickname.
    
    Returns (wavelength, transmission) vectors.
    """

    filter_directory = load_filters()
    filter_nicknames = Dict{String, String}()
    filter_names = sort(collect(keys(filter_directory)))

    for key in filter_names
        for n in filter_directory[key]["nicknames"]
            filter_nicknames[n] = key
        end
    end

    if haskey(filter_nicknames, nickname)
        real_name = filter_nicknames[nickname]
        filt = readdlm(filterpath * real_name)
        return filt[:,1], filt[:,2]
    else
        throw(LazyError("Filter '$nickname' not found"))
    end
end

function get_pivot_wavelengths(bands::Vector{String})::Vector{Float64}
    """
    Compute the pivot wavelengths for the given bands.
    
    Returns vector of pivot wavelengths in Angstroms.
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
    template_directory = load_templates()
    template_set_names = sort(collect(keys(template_directory)))
    for template_set_name in template_set_names
        template_files = template_directory[template_set_name]["files"]
        for template_file in template_files
            if !isfile(templatepath * template_file)
                throw(LazyError("Template file $template_file not found"))
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

function params()
    """ 
    Read in the example parameter file and print the parameters.
    """
    param_file = @__DIR__() * "/example_params.toml"
    if !isfile(param_file)
        throw(LazyError("Parameter file $param_file not found"))
    end
    f = readlines(param_file)
    for line in f
        println(line)
    end


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
            throw(LazyError("old_fluxes and old_errs must be the same length"))
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


function load_template(file::String)::Tuple{Vector{Float64}, Union{Vector{Float64}, Matrix{Float64}}, Union{Nothing, Vector{Float64}}}
    """
    Load a template from a file.
    
    Returns (wavelength, flux, redshift_grid) where redshift_grid is Nothing for non-redshift-dependent templates.
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
        else
            templz = nothing
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

