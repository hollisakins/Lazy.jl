"""
Photometric redshift fitting functions for Lazy.jl

Contains the main fit() function, fit_streaming() implementation, single object fitting,
and template error precomputation.
"""

using TOML
using FITSIO
using NonNegLeastSquares
using LinearAlgebra
using Base.Threads
using Trapz

# Flux zero points for AB magnitude system: m_AB = -2.5*log10(F_cgs) - 48.6
const FLUX_ZP = Dict(
    "cgs" => -48.6, "Jy" => 8.9,
    "uJy" => 23.9, "ujy" => 23.9,
    "nJy" => 31.4, "njy" => 31.4,
)

const RF_MAG_NAMES = ["M_UV", "M_U", "M_V", "M_J"]

"""
    luminosity_distance_Mpc(z, H0, Om)

Luminosity distance in Mpc for flat LCDM cosmology using Simpson's rule.
"""
function luminosity_distance_Mpc(z::Float64, H0::Float64, Om::Float64)::Float64
    z <= 0.0 && return 0.0
    c_km_s = 299792.458
    OL = 1.0 - Om
    n = 500
    dz = z / n
    integral = 0.0
    for k in 0:n
        zp = k * dz
        Ezp = sqrt(Om * (1.0 + zp)^3 + OL)
        w = (k == 0 || k == n) ? 1.0 : (k % 2 == 1 ? 4.0 : 2.0)
        integral += w / Ezp
    end
    integral *= dz / 3.0
    return (1.0 + z) * (c_km_s / H0) * integral
end

"""
    distance_modulus(z, H0, Om)

Distance modulus DM = 5*log10(d_L_Mpc) + 25.
"""
function distance_modulus(z::Float64, H0::Float64, Om::Float64)::Float64
    d_L = luminosity_distance_Mpc(z, H0, Om)
    d_L <= 0.0 && return 0.0
    return 5.0 * log10(d_L) + 25.0
end

function fit_single_object(j::Int, fnu_j::Vector{Float64}, efnu_j::Vector{Float64},
                          templgrid::Array{Float64,3}, template_error_grid::Matrix{Float64},
                          zgrid::Vector{Float64}, nphot_min::Int, nband::Int, ntempl::Int, nz::Int;
                          interrupted_flag::Union{Nothing,Ref{Bool}}=nothing,
                          z_fix_idx::Int=-1)
    
    # Wrap entire function in try-catch to handle individual object failures gracefully
    try
        # Pre-allocate working arrays to avoid repeated allocations
        templgrid_ij = Matrix{Float64}(undef, nband, ntempl)
        fnu_mod_j = Vector{Float64}(undef, nband)
        chi2_row = Vector{Float32}(undef, nz)

        # Handle fixed redshift (z_spec mode)
        if z_fix_idx > 0
            fill!(chi2_row, Float32(1e18))
            z_range = z_fix_idx:z_fix_idx
        else
            z_range = 1:nz
        end

        # Track best fit for this object
        best_chi2 = Inf
        best_z_idx = 0
        best_coeffs = zeros(ntempl)

        # Loop over redshifts for this object
        for i in z_range
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
        chi2_row_failed = fill(Float32(-1.0), nz)
        
        return (zbest_failed, chi2best_failed, coeffsbest_failed, chi2_row_failed)
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

"""
    finalize_output(work_file, output_file, output_format, preserve_work_file)

Handle final output dispatch based on output_format ("fits", "hdf5", or "both").
"""
function finalize_output(work_file::String, output_file::String, output_format::String, preserve_work_file::Bool)
    if output_format == "fits" || output_format == "both"
        # Derive FITS filename from output_file
        fits_file = endswith(output_file, ".fits") ? output_file : splitext(output_file)[1] * ".fits"
        println("💾 Exporting to FITS: $fits_file")
        convert_hdf5_to_fits(work_file, fits_file)
    end

    if output_format == "hdf5" || output_format == "both"
        # Derive HDF5 filename from output_file
        hdf5_file = if endswith(output_file, ".h5") || endswith(output_file, ".hdf5")
            output_file
        else
            splitext(output_file)[1] * ".h5"
        end

        if output_format == "both"
            # Copy work file (don't move — still need it or preserve_work_file may keep it)
            cp(work_file, hdf5_file, force=true)
        else
            mv(work_file, hdf5_file, force=true)
        end
        println("✅ Results saved to: $hdf5_file")
    end

    # Clean up work file if it still exists and we're not preserving it
    if isfile(work_file)
        if preserve_work_file
            println("💾 Work file preserved: $work_file")
        else
            rm(work_file)
            println("🗑️ Work file automatically removed: $work_file")
        end
    end
end

function fit(param)

    # Load in TOML parameter file
    if !isfile(param)
        throw(LazyError("parameter file $param not found"))
    end
    if !endswith(param, ".toml")
        throw(LazyError("parameter file $param must have a .toml extension"))
    end
    param = with_spinner(() -> TOML.parsefile(param), "Loading parameter file: $param", "📊")
    

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

        # Rest-frame absolute magnitudes
        output_restframe_mags = get(io, "output_restframe_mags", false)
        flux_units = nothing
        flux_zp = 0.0
        cosmo_H0 = 70.0
        cosmo_Om = 0.3
        if output_restframe_mags
            if !haskey(io, "flux_units")
                throw(LazyError("parameter `io.flux_units` is required when output_restframe_mags = true"))
            end
            flux_units = io["flux_units"]
            if !haskey(FLUX_ZP, flux_units)
                throw(LazyError("Unrecognized flux_units '$flux_units'. Must be one of: cgs, Jy, uJy, nJy"))
            end
            flux_zp = FLUX_ZP[flux_units]
            cosmo_H0 = Float64(get(io, "H0", 70.0))
            cosmo_Om = Float64(get(io, "Om", 0.3))
            println("📊 Rest-frame magnitudes: enabled (flux_units=$flux_units, H0=$cosmo_H0, Om=$cosmo_Om)")
        end

        # Parse output format: "fits", "hdf5", or "both"
        # If not specified, infer from output_file extension
        output_format = get(io, "output_format", nothing)
        if output_format === nothing
            if endswith(output_file, ".h5") || endswith(output_file, ".hdf5")
                output_format = "hdf5"
            else
                output_format = "fits"
            end
        end
        output_format = lowercase(output_format)
        if !(output_format in ["fits", "hdf5", "both"])
            throw(LazyError("parameter `io.output_format` must be 'fits', 'hdf5', or 'both' (got '$output_format')"))
        end
    else
        throw(LazyError("section `io` not found in the parameter file"))
    end

    # Parse runtime configuration
    chunked_processing = false
    target_memory_gb = 0.5
    preserve_work_file = false
    
    if haskey(param, "runtime")
        runtime = param["runtime"]
        
        if haskey(runtime, "chunked_processing")
            chunked_processing = runtime["chunked_processing"]
        end
        
        if haskey(runtime, "target_memory_gb")
            target_memory_gb = runtime["target_memory_gb"]
            if target_memory_gb <= 0
                throw(LazyError("runtime.target_memory_gb must be positive"))
            end
        end
        
        if haskey(runtime, "preserve_work_file")
            preserve_work_file = runtime["preserve_work_file"]
        end
        
        # # Provide user feedback about runtime configuration
        # if chunked_processing
        #     println("⚙️ Runtime: Chunked processing enabled")
        #     println("   Target memory: $(target_memory_gb) GB per chunk")
        #     if preserve_work_file
        #         println("   Work file preservation: enabled")
        #     end
        # else
        #     println("⚙️ Runtime: In-memory processing (default)")
        # end
        # else
        # println("⚙️ Runtime: In-memory processing (default)")
    end

    IDs = Int.(read(cat[2], "ID"))
    nobj = length(IDs)
    println("📊 Objects: $nobj")
        
    println("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
    if !haskey(param, "translate")
        throw(LazyError("section `translate` not found in the parameter file"))
    end
    
    translate = param["translate"]
    # Extract z_spec column mapping before building bands list
    zspec_column = haskey(translate, "zspec") ? pop!(translate, "zspec") : nothing
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

        


    println("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")

    
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

    println("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")

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

    # Estimate and report memory usage
    memory_estimate = estimate_memory_usage(nobj, nz, nband, ntempl)
    print_memory_estimate(memory_estimate; chunked_processing=chunked_processing, target_memory_gb=target_memory_gb)
    
    # Check for template cache before building grid
    use_cache = get(fitting, "template_cache", true)  # Default enabled

    # CGM Lyman-alpha damping wing model (Asada+24)
    add_cgm = get(fitting, "add_cgm", true)
    cgm_A = Float64(get(fitting, "cgm_A", 3.5918))
    cgm_a = Float64(get(fitting, "cgm_a", 1.8414))
    cgm_c = Float64(get(fitting, "cgm_c", 18.001))
    if add_cgm
        println("🌌 CGM damping wing: enabled (Asada+24)")
    else
        println("🌌 CGM damping wing: disabled")
    end

    # Spectroscopic redshift fixing
    use_zspec = get(fitting, "use_zspec", false)
    zspec = fill(NaN, nobj)
    if use_zspec
        if zspec_column !== nothing
            all_cols = FITSIO.colnames(cat[2])
            if zspec_column in all_cols
                zspec = Float64.(read(cat[2], zspec_column))
                zspec[.!isfinite.(zspec) .| (zspec .< 0)] .= NaN
                nzspec = sum(isfinite.(zspec))
                println("🔭 Spectroscopic redshift fixing: enabled ($nzspec / $nobj objects with z_spec)")
            else
                @warn "use_zspec=true but column '$zspec_column' not found in catalog. Fitting all objects with full grid."
                use_zspec = false
            end
        else
            @warn "use_zspec=true but no 'zspec' entry in [translate] section. Fitting all objects with full grid."
            use_zspec = false
        end
    end

    templgrid = nothing
    template_error_grid = nothing
    cache_key = ""
    
    if use_cache
        # Generate cache key from all relevant parameters
        cache_key = generate_cache_key(templates, zgrid, bands, fitting["igm_model"],
                                     fitting["template_error"], fitting["template_error_scale"];
                                     add_cgm=add_cgm, cgm_A=cgm_A, cgm_a=cgm_a, cgm_c=cgm_c)
        
        # Try to load from cache
        cached_data = load_template_cache(cache_key)
        if cached_data !== nothing
            templgrid, template_error_grid = cached_data
            println("✅ Loading cached template grid ($(cache_key)) - " * summary(templgrid))
        end
    end
    
    # Build template grid if not loaded from cache
    if templgrid === nothing
        if !haskey(fitting, "igm_model")
            throw(LazyError("parameter `igm_model` not found in the parameter file"))
        end
        
        # Use unified template grid builder for photometry mode
        templgrid, template_error_grid_new, _ = build_template_grid(
            templates, zgrid, fitting["igm_model"], fitting;
            output_type=:photometry, bands=bands,
            add_cgm=add_cgm, cgm_A=cgm_A, cgm_a=cgm_a, cgm_c=cgm_c
        )
        
        # Use newly built template error grid if not loaded from cache
        if template_error_grid === nothing  
            template_error_grid = template_error_grid_new
        end
        
        # Report template grid size
        templgrid_size_gb = sizeof(templgrid) / (1024^3)
        println("✅ Template grid built: " * summary(templgrid) * " (~$(round(templgrid_size_gb, digits=2)) GB)")
        if templgrid_size_gb > 8.0
            println("⚠️ Large template grid detected. Consider reducing nz, nband, or ntempl if memory issues occur.")
        end
        
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
                "nband" => nband,
                "add_cgm" => add_cgm,
                "cgm_A" => cgm_A,
                "cgm_a" => cgm_a,
                "cgm_c" => cgm_c
            )
            
            save_template_cache(templgrid, template_error_grid, metadata, cache_key)
            println("💾 Template grid cached for future use ($(cache_key))")
        end
    end

    println("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")

    # Template grids are now ready (either from cache or newly built)

    # Build rest-frame template grid if needed
    restframe_templgrid = nothing
    if output_restframe_mags
        restframe_templgrid = with_spinner(
            () -> build_restframe_template_grid(templates, zgrid),
            "Building rest-frame template grid", "📊")
        println("✅ Rest-frame template grid built: " * summary(restframe_templgrid))
    end

    # Always use fit_streaming() - it handles both chunked and in-memory modes
    return fit_streaming(param, templgrid, template_error_grid, zgrid, bands,
                       templates, nobj, nz, nband, ntempl, fnu, efnu, IDs,
                       nphot_min, output_file, output_pz, output_templ,
                       target_memory_gb, preserve_work_file, chunked_processing,
                       use_zspec, zspec, output_format;
                       output_restframe_mags=output_restframe_mags,
                       restframe_templgrid=restframe_templgrid,
                       flux_zp=flux_zp, cosmo_H0=cosmo_H0, cosmo_Om=cosmo_Om)
end

"""
    fit_streaming(param::Dict, templgrid::Array{Float64,3}, template_error_grid::Matrix{Float64},
                  zgrid::Vector{Float64}, bands::Vector{String}, templates::Vector{String},
                  nobj::Int, nz::Int, nband::Int, ntempl::Int,
                  fnu::Matrix{Float64}, efnu::Matrix{Float64}, IDs::Vector{Int},
                  nphot_min::Int, output_file::String, output_pz::Bool, output_templ::Bool,
                  target_memory_gb::Float64, preserve_work_file::Bool, chunked_processing::Bool,
                  use_zspec::Bool, zspec::Vector{Float64}, output_format::String="fits")

Main streaming fit function that handles both chunked and in-memory processing modes.
When chunked_processing=false, uses single chunk (in-memory mode).
When chunked_processing=true, uses memory-controlled chunking.
output_format controls final output: "fits", "hdf5", or "both".
"""
function fit_streaming(param::Dict, templgrid::Array{Float64,3}, template_error_grid::Matrix{Float64},
                      zgrid::Vector{Float64}, bands::Vector{String}, templates::Vector{String},
                      nobj::Int, nz::Int, nband::Int, ntempl::Int,
                      fnu::Matrix{Float64}, efnu::Matrix{Float64}, IDs::Vector{Int},
                      nphot_min::Int, output_file::String, output_pz::Bool, output_templ::Bool,
                      target_memory_gb::Float64, preserve_work_file::Bool, chunked_processing::Bool,
                      use_zspec::Bool, zspec::Vector{Float64}, output_format::String="fits";
                      output_restframe_mags::Bool=false,
                      restframe_templgrid::Union{Nothing,Array{Float64,3}}=nothing,
                      flux_zp::Float64=0.0, cosmo_H0::Float64=70.0, cosmo_Om::Float64=0.3)
    
    # CGM params (read from fitting section)
    fitting = param["fitting"]
    add_cgm = get(fitting, "add_cgm", true)
    cgm_A = Float64(get(fitting, "cgm_A", 3.5918))
    cgm_a = Float64(get(fitting, "cgm_a", 1.8414))
    cgm_c = Float64(get(fitting, "cgm_c", 18.001))

    # Generate work file name
    work_file = output_file * ".work.h5"
    
    # Check for resume capability
    can_resume, last_obj = check_resume_file(work_file)
    start_obj = 1
    
    if can_resume
        action, should_resume = prompt_resume(work_file, last_obj, nobj)
        
        if action == "resume" && should_resume
            start_obj = last_obj + 1
            println("✅ Resuming from object $start_obj")
        elseif action == "complete"
            println("✅ Marking run as complete and generating output...")
            finalize_hdf5_work_file(work_file)
            
            # Jump to output conversion
            finalize_output(work_file, output_file, output_format, preserve_work_file)
            return 0  # Exit early for complete runs
        elseif action == "keep"
            println("💾 Keeping work file and exiting...")
            return 0  # Exit without processing
        else  # action == "overwrite"
            println("🔄 Starting fresh run")
            rm(work_file, force=true)
            can_resume = false
        end
    end
    
    # Create work file if not resuming
    if !can_resume
        println("📝 Creating HDF5 work file: $work_file")
        create_hdf5_work_file(work_file, param, nobj, nz, nband, ntempl, bands, zgrid, templates)
    end

    # Write z_spec data to work file if available
    if use_zspec
        h5open(work_file, "r+") do file
            if haskey(file["results"], "z_spec")
                file["results/z_spec"][:] = zspec
            end
        end
    end
    
    # Determine chunk size based on processing mode
    if chunked_processing
        # Use memory-controlled chunking
        chunk_size = get_chunk_size_for_memory_target(nobj, nz, target_memory_gb)
        actual_memory_gb = chunk_size * nz * 4 / (1024^3)  # Float32 chi2 values
        println("⚙️ Chunked processing: $(format_number(chunk_size)) objects per chunk (~$(round(actual_memory_gb, digits=2)) GB per chunk)")
    else
        # In-memory mode: single chunk = entire catalog
        chunk_size = nobj
        estimated_gb = nobj * nz * 4 / (1024^3)  # Rough estimate for display
        println("⚙️ In-memory processing: all $(format_number(nobj)) objects in single batch (~$(round(estimated_gb, digits=2)) GB)")
        if estimated_gb > 8.0
            println("   ⚠️  High memory usage detected. Consider setting chunked_processing = true if you encounter memory issues")
        end
    end
    
    # Initialize progress tracking
    total_remaining = nobj - start_obj + 1
    progress_bar = ProgressBar(total_remaining, "Fitting objects...", 
                              show_rate=true, show_eta=true, prefix_emoji="🧠")
    
    # Display initial progress at 0%
    display_progress(progress_bar)
    
    # Process objects in chunks
    objects_processed = 0
    
    # Set up interrupt handling
    interrupted = Ref(false)
    
    try
        for chunk_start in start_obj:chunk_size:nobj
            chunk_end = min(chunk_start + chunk_size - 1, nobj)
            chunk_nobj = chunk_end - chunk_start + 1
            
            # Process this chunk
            chunk_IDs = IDs[chunk_start:chunk_end]
            chunk_fnu = fnu[chunk_start:chunk_end, :]
            chunk_efnu = efnu[chunk_start:chunk_end, :]
            
            # Initialize chunk result arrays
            chunk_zbest = zeros(chunk_nobj)
            chunk_chi2best = zeros(chunk_nobj)
            chunk_coeffsbest = zeros(chunk_nobj, ntempl)
            chunk_chi2grid = zeros(Float32, chunk_nobj, nz)
            
            # Fit objects in this chunk using threading
            nthreads_to_use = min(Threads.nthreads(), chunk_nobj)
            objects_per_task = max(1, chunk_nobj ÷ nthreads_to_use)
            ntasks = cld(chunk_nobj, objects_per_task)
            
            tasks = Vector{Task}(undef, ntasks)
            for task_id in 1:ntasks
                task_start = (task_id - 1) * objects_per_task + 1
                task_end = min(task_id * objects_per_task, chunk_nobj)
                
                tasks[task_id] = Threads.@spawn begin
                    batch_results = Vector{Tuple{Int, Float64, Float64, Vector{Float64}, Vector{Float32}}}()
                    
                    for j_local in task_start:task_end
                        # Check for interruption before starting work
                        if interrupted[]
                            push!(batch_results, (j_local, -1.0, -1.0, zeros(ntempl), fill(Float32(-1.0), nz)))
                            increment!(progress_bar)
                            continue
                        end
                        
                        j_global = chunk_start + j_local - 1

                        # Compute fixed redshift index if z_spec available
                        z_fix_idx_j = -1
                        if use_zspec && isfinite(zspec[j_global])
                            z_fix_idx_j = argmin(abs.(zgrid .- zspec[j_global]))
                        end

                        zbest_j, chi2best_j, coeffsbest_j, chi2_row_j = fit_single_object(
                            j_global, chunk_fnu[j_local,:], chunk_efnu[j_local,:],
                            templgrid, template_error_grid, zgrid, nphot_min,
                            nband, ntempl, nz;
                            interrupted_flag=interrupted, z_fix_idx=z_fix_idx_j
                        )
                        
                        push!(batch_results, (j_local, zbest_j, chi2best_j, coeffsbest_j, chi2_row_j))
                        increment!(progress_bar)
                    end
                    
                    return batch_results
                end
            end
            
            # Collect results from tasks
            for task in tasks
                batch_results = fetch(task)
                for (j_local, zbest_j, chi2best_j, coeffsbest_j, chi2_row_j) in batch_results
                    chunk_zbest[j_local] = zbest_j
                    chunk_chi2best[j_local] = chi2best_j
                    chunk_coeffsbest[j_local, :] = coeffsbest_j
                    chunk_chi2grid[j_local, :] = chi2_row_j
                end
            end
            
            # Calculate P(z) statistics for this chunk (always computed)
            pz_chunk = exp.(-0.5 * chunk_chi2grid)
            cpz_chunk = cumsum(pz_chunk, dims=2) ./ sum(pz_chunk, dims=2)

            chunk_z_l95 = zgrid[map(argmin, eachrow(abs.(cpz_chunk .- 0.025)))]
            chunk_z_l68 = zgrid[map(argmin, eachrow(abs.(cpz_chunk .- 0.160)))]
            chunk_z_med = zgrid[map(argmin, eachrow(abs.(cpz_chunk .- 0.500)))]
            chunk_z_u68 = zgrid[map(argmin, eachrow(abs.(cpz_chunk .- 0.840)))]
            chunk_z_u95 = zgrid[map(argmin, eachrow(abs.(cpz_chunk .- 0.975)))]

            # Integrated P(z) quantities: P(z > Z) and Sz
            z_integers = collect(1:floor(Int, maximum(zgrid)))
            chunk_pz_gt = fill(Float32(-1), chunk_nobj, length(z_integers))
            z_max_grid = maximum(zgrid)
            Sz_bc = collect(0.0:1.0:z_max_grid)
            chunk_Sz = fill(-1.0, chunk_nobj)

            for j in 1:chunk_nobj
                total = trapz(zgrid, @view pz_chunk[j, :])
                total <= 0 && continue

                # P(z > Z) for each integer Z
                for (zi, Z) in enumerate(z_integers)
                    idx = searchsortedfirst(zgrid, Z + eps())
                    idx > length(zgrid) && continue
                    chunk_pz_gt[j, zi] = Float32(trapz(@view(zgrid[idx:end]), @view(pz_chunk[j, idx:end])) / total)
                end

                # Sz: center of the unit-width bin with the most P(z)
                best_prob = -1.0
                best_bc = -1.0
                for bc in Sz_bc
                    zlo = bc == 0.0 ? 0.0 : bc - 0.5
                    zhi = bc == z_max_grid ? z_max_grid : bc + 0.5
                    ilo = searchsortedfirst(zgrid, zlo)
                    ihi = searchsortedlast(zgrid, zhi - eps())
                    ilo > ihi && continue
                    prob = trapz(@view(zgrid[ilo:ihi]), @view(pz_chunk[j, ilo:ihi])) / total
                    if prob > best_prob
                        best_prob = prob
                        best_bc = bc
                    end
                end
                chunk_Sz[j] = best_bc
            end

            # Only save full chi2 grid when output_pz is requested
            if !output_pz
                chunk_chi2grid = nothing
            end
            
            # Calculate best-fit photometry for this chunk
            chunk_photobest = zeros(chunk_nobj, nband)
            for j_local in 1:chunk_nobj
                if chunk_zbest[j_local] > 0  # Only for valid fits
                    z_idx = argmin(abs.(zgrid .- chunk_zbest[j_local]))
                    for k in 1:nband
                        chunk_photobest[j_local, k] = 0.0
                        for t in 1:ntempl
                            chunk_photobest[j_local, k] += templgrid[t, z_idx, k] * chunk_coeffsbest[j_local, t]
                        end
                    end
                end
            end
            
            # Compute rest-frame absolute magnitudes for this chunk
            chunk_restframe_mags = nothing
            if output_restframe_mags && restframe_templgrid !== nothing
                nrf = 4
                chunk_restframe_mags = fill(NaN, chunk_nobj, nrf)
                for j_local in 1:chunk_nobj
                    chunk_zbest[j_local] <= 0 && continue
                    z_idx = argmin(abs.(zgrid .- chunk_zbest[j_local]))
                    DM = distance_modulus(chunk_zbest[j_local], cosmo_H0, cosmo_Om)
                    kcorr = 2.5 * log10(1.0 + chunk_zbest[j_local])
                    for k in 1:nrf
                        f_rest = 0.0
                        for t in 1:ntempl
                            f_rest += restframe_templgrid[t, z_idx, k] * chunk_coeffsbest[j_local, t]
                        end
                        if f_rest > 0.0
                            chunk_restframe_mags[j_local, k] = -2.5 * log10(f_rest) + flux_zp - DM + kcorr
                        end
                    end
                end
            end

            # Write chunk results to HDF5
            write_chunk_results(work_file, chunk_start, chunk_end,
                              chunk_IDs, chunk_zbest, chunk_chi2best, chunk_coeffsbest,
                              chunk_z_l95, chunk_z_l68, chunk_z_med, chunk_z_u68, chunk_z_u95,
                              chunk_photobest, bands, chunk_chi2grid, zgrid;
                              pz_gt=chunk_pz_gt, Sz=chunk_Sz, z_integers=z_integers,
                              restframe_mags=chunk_restframe_mags)
            
            # Progress is now updated per-object within the threaded tasks
            objects_processed += chunk_nobj
        end
        
        # Finish progress bar immediately after fitting loop completes
        if interrupted[]
            finish!(progress_bar, final_message="interrupted")
        else
            finish!(progress_bar)
        end
        
        # Generate template output if requested (after fitting is complete)
        if output_templ
            println("📊 Generating template output...")
            
            # Use unified template grid builder for spectral mode
            templgrid_spectral, _, wavelength_grid = build_template_grid(
                templates, zgrid, param["fitting"]["igm_model"], param["fitting"];
                output_type=:spectral,
                add_cgm=add_cgm, cgm_A=cgm_A, cgm_a=cgm_a, cgm_c=cgm_c
            )
            
            # Write template data to HDF5 work file
            write_template_data(work_file, templgrid_spectral, wavelength_grid, zgrid, templates)
            println("✅ Template output generated")
        end
        
        # Mark work file as complete
        finalize_hdf5_work_file(work_file)
        
        # Handle output format and work file preservation
        finalize_output(work_file, output_file, output_format, preserve_work_file)
        
    catch e
        if isa(e, InterruptException)
            interrupted[] = true
            println("\n⚠️ Fitting interrupted by user. Processing results for completed objects...")
            finalize_hdf5_work_file(work_file)  # Mark as complete for resume
            finish!(progress_bar, final_message="interrupted")
        else
            println("\n⚠️ Error during streaming fit: $e")
            println("💾 Partial results preserved in: $work_file")
            rethrow(e)
        end
    end
    
    return 0
end