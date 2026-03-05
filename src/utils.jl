"""
Utility functions for Lazy.jl

Contains progress bars, spinners, formatting utilities, memory estimation, and error handling.
"""

using Base.Threads

# Custom exception type for Lazy.jl
struct LazyError <: Exception
    msg::String
    LazyError(msg::String) = new(msg)
end

Base.showerror(io::IO, e::LazyError) = print(io, "LazyError: ", e.msg)

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
                        bar_width::Int=30, prefix_emoji::String="🧠")
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
        elseif current == 0
            # Show initial message when no progress yet
            if pb.show_rate || pb.show_eta
                display_str *= "; starting..."
            end
        end
    end
    
    display_str *= ")    " # Ensure trailing spaces to avoid duplicate characters
    
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
        
        # Run operation on a spawned thread, keep spinner on main thread
        result_ref = Ref{Any}(nothing)
        error_ref = Ref{Any}(nothing)
        done = Ref(false)

        work_task = Threads.@spawn begin
            try
                result_ref[] = operation()
            catch e
                error_ref[] = e
            finally
                done[] = true
            end
        end

        # Spinner runs on main thread — guaranteed terminal access
        i = 1
        while !done[]
            spinner_line = "$(SPINNER_CHARS[i]) $message..."
            print("\r$spinner_line")
            max_line_length[] = max(max_line_length[], length(spinner_line))
            flush(stdout)
            sleep(0.1)
            i = (i % length(SPINNER_CHARS)) + 1
        end
        wait(work_task)

        # Clear spinner and print result
        clear_line = "\r" * " "^max_line_length[] * "\r"
        if error_ref[] !== nothing
            print(clear_line * "❌ $message (failed)\n")
            flush(stdout)
            throw(error_ref[])
        else
            print("$clear_line$success_emoji $message\n")
            flush(stdout)
            return result_ref[]
        end
    else
        # Fallback for non-interactive terminals - only print completion message
        result = operation()
        println("$success_emoji $message")
        return result
    end
end

function estimate_memory_usage(nobj::Int, nz::Int, nband::Int, ntempl::Int)::Dict{String, Float64}
    """
    Estimate memory usage for the fitting process in GB.
    
    Returns a dictionary with estimated memory usage for different components.
    """
    # Template grid: nband × ntempl × nz × 4 bytes (Float32)
    templgrid_gb = (ntempl * nz * nband * 4) / (1024^3)
    
    # Chi2 grid: nobj × nz × 4 bytes (Float32)  
    chi2grid_gb = (nobj * nz * 4) / (1024^3)
    
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

function print_memory_estimate(mem_dict::Dict{String, Float64}; chunked_processing::Bool=false, target_memory_gb::Float64=0.5)
    """
    Print a formatted memory usage estimate.
    """
    if chunked_processing
        peak_gb = mem_dict["estimated_peak"]
        println("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
        println("📊 Memory: ~$(round(peak_gb, digits=2)) GB peak (chunked processing enabled, ~$(target_memory_gb) GB per chunk)")
        
        if peak_gb > 32.0
            println("   ⚠️  WARNING: Very high memory usage detected! Consider reducing chunk size")
        elseif peak_gb > 16.0  
            println("   ⚠️  CAUTION: High memory usage detected!")
        end
        println("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
    else
        peak_gb = mem_dict["estimated_peak"]
        println("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
        println("📊 Memory: ~$(round(peak_gb, digits=2)) GB peak")
        
        if peak_gb > 32.0
            println("   ⚠️  WARNING: Very high memory usage detected! Consider using chunked processing")
        elseif peak_gb > 16.0  
            println("   ⚠️  CAUTION: High memory usage detected! Consider using chunked processing")
        end
        println("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
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