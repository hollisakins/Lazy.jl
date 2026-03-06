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

function with_spinner(operation::Function, message::String, success_emoji::String="✓")
    """
    Run an operation with a spinner for visual feedback.
    Falls back to static messages in non-interactive terminals.
    """
    t0 = time()
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
        elapsed = format_time(time() - t0)
        clear_line = "\r" * " "^max_line_length[] * "\r"
        if error_ref[] !== nothing
            print(clear_line * "❌ $message (failed)\n")
            flush(stdout)
            throw(error_ref[])
        else
            print("$clear_line$success_emoji $message ($elapsed)\n")
            flush(stdout)
            return result_ref[]
        end
    else
        # Fallback for non-interactive terminals - only print completion message
        result = operation()
        elapsed = format_time(time() - t0)
        println("$success_emoji $message ($elapsed)")
        return result
    end
end

function estimate_peak_memory(nobj::Int, nz::Int, nband::Int, ntempl::Int, chunk_size::Int,
                              zmax::Float64;
                              output_pz::Bool=false, output_restframe_mags::Bool=false,
                              output_forced_lowz::Bool=false)
    """
    Estimate peak memory usage in GB.

    Persistent arrays (live for entire run) + per-chunk arrays (one chunk at a time).
    Peak occurs during P(z) computation when fitting results and P(z) stats coexist.
    """
    n_z_integers = max(1, floor(Int, zmax))

    # === Persistent arrays (always in memory) ===
    # Template grid: nband × ntempl × nz × Float32
    templgrid_bytes = nband * ntempl * nz * 4
    # Template error grid: nz × nband × Float64
    template_error_bytes = nz * nband * 8
    # Input photometry: fnu + efnu (nobj × nband × Float64 × 2)
    photdata_bytes = nobj * nband * 2 * 8
    # IDs + zspec: nobj × 8 × 2
    ids_bytes = nobj * 8 * 2
    # Rest-frame template grid: 4 × ntempl × nz × Float64 (if enabled)
    restframe_bytes = output_restframe_mags ? 4 * ntempl * nz * 8 : 0

    persistent_bytes = templgrid_bytes + template_error_bytes + photdata_bytes + ids_bytes + restframe_bytes

    # === Per-chunk arrays (one chunk live at a time) ===
    C = chunk_size
    # Fitting results: chi2grid (dominant), zbest, chi2best, coeffs, lowz variants
    chi2grid_bytes = nz * C * 4  # Float32
    fitting_scalars_bytes = C * 8 * 2  # zbest + chi2best
    coeffs_bytes = C * ntempl * 8 * 2  # coeffsbest + coeffsbest_lowz
    lowz_scalars_bytes = C * 8 * 2  # zbest_lowz + chi2best_lowz

    # P(z) statistics arrays (coexist with fitting results at peak)
    pz_quantiles_bytes = C * 8 * 5  # z_l95, z_l68, z_med, z_u68, z_u95
    pz_gt_bytes = C * n_z_integers * 4  # Float32
    pz_bins_bytes = C * (n_z_integers + 1) * 4  # n_bins ≈ n_z_integers + 1
    pz_misc_bytes = C * (8 + 4 + 4)  # Sz(Float64) + pz_cen(Float32) + pz_zgtrzb2(Float32)
    pz_grid_bytes = output_pz ? nz * C * 4 : 0  # Float32, same size as chi2grid

    # Forced low-z P(z) arrays
    lowz_pz_bytes = output_forced_lowz ? C * (8 * 6 + nband * 8) : 0  # quantiles + delta_chi2 + photobest_lowz

    # Best-fit photometry
    photobest_bytes = C * nband * 8
    restframe_mags_bytes = output_restframe_mags ? C * 4 * 8 : 0

    chunk_bytes = chi2grid_bytes + fitting_scalars_bytes + coeffs_bytes + lowz_scalars_bytes +
                  pz_quantiles_bytes + pz_gt_bytes + pz_bins_bytes + pz_misc_bytes + pz_grid_bytes +
                  lowz_pz_bytes + photobest_bytes + restframe_mags_bytes

    peak_bytes = persistent_bytes + chunk_bytes
    peak_gb = peak_bytes / (1024^3)

    return Dict(
        "persistent" => persistent_bytes / (1024^3),
        "per_chunk" => chunk_bytes / (1024^3),
        "estimated_peak" => peak_gb
    )
end

function print_memory_estimate(peak_gb::Float64; chunked::Bool=false)
    println("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
    println("📊 Estimated peak memory: ~$(round(peak_gb, digits=2)) GB")
    if peak_gb > 32.0
        msg = chunked ? "Consider reducing chunk size" : "Consider using chunked processing"
        println("   ⚠️  WARNING: Very high memory usage detected! $msg")
    elseif peak_gb > 16.0
        msg = chunked ? "" : " Consider using chunked processing"
        println("   ⚠️  CAUTION: High memory usage detected!$msg")
    end
    println("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
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