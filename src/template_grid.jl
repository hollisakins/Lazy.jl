"""
Template Grid Generation for Lazy.jl

Unified template grid builder supporting both photometry and spectral output modes.
Eliminates code duplication between photometric fitting and template export.
"""

using HDF5
using LinearAlgebra
using Interpolations
using Base.Threads
using DelimitedFiles
using Trapz

"""
    determine_common_wavelength_grid(templates::Vector{String})

Find the common wavelength grid by selecting the template with the largest wavelength array.
This is used for spectral template output to ensure all templates are on the same grid.

Returns: wavelength grid from the highest resolution template
"""
function determine_common_wavelength_grid(templates::Vector{String})
    common_templwav = nothing
    
    for template in templates
        templwav_i, templfnu_i, templz_i = load_template(template)
        if common_templwav === nothing
            common_templwav = templwav_i
        elseif length(templwav_i) > length(common_templwav)
            common_templwav = templwav_i
        end
    end
    
    return common_templwav
end

"""
    load_igm_model(igm_model_name::String)

Load IGM transmission data from HDF5 file.

Returns: (igm_redshifts, igm_wavelengths, igm_transmission, idx_igm, maxz)
"""
function load_igm_model(igm_model_name::String)
    igm_file = h5open(igmpath * igm_model_name * ".hdf5", "r")
    igm_redshifts = igm_file["redshifts"][:]
    igm_wavelengths = igm_file["wavelengths"][:]
    igm_transmission = igm_file["transmission"][:,:]
    close(igm_file)
    
    idx_igm = searchsortedfirst(igm_wavelengths, 1215.67)
    maxz = igm_redshifts[end]
    
    return igm_redshifts, igm_wavelengths, igm_transmission, idx_igm, maxz
end

"""
    sigma_a(nu_rest)
Ly-alpha absorption cross-section (Lorentzian damping wing). Returns cm².
"""
function sigma_a(nu_rest::Float64)::Float64
    Lam_a = 6.255486e8
    nu_lya = 2.46607e15
    C = 6.9029528e22
    nu_ratio = nu_rest / nu_lya
    sig = C * nu_ratio^4 / (4.0 * π^2 * (nu_rest - nu_lya)^2 + Lam_a^2 * nu_ratio^6 / 4.0)
    return sig * 1e-16
end

"""
    cgm_sigmoid(z, A, a, c)
HI column density vs redshift sigmoid (Asada+24). Returns log10(N_HI).
"""
function cgm_sigmoid(z::Float64, A::Float64, a::Float64, c::Float64)::Float64
    return A / (1.0 + exp(-a * (z - 6.0))) + c
end

"""
    apply_cgm_damping_wing!(transmission, templwav, z, cgm_A, cgm_a, cgm_c)
In-place multiply transmission by exp(-τ_CGM). No-op for z < 6.
"""
function apply_cgm_damping_wing!(transmission::Vector{Float64}, templwav::Vector{Float64},
                                  z::Float64, cgm_A::Float64, cgm_a::Float64, cgm_c::Float64)
    z < 6.0 && return
    N_HI = 10.0^cgm_sigmoid(z, cgm_A, cgm_a, cgm_c)
    c_ang = 2.99792458e18
    @inbounds for i in eachindex(transmission)
        nu_rest = c_ang / templwav[i]
        transmission[i] *= exp(-N_HI * sigma_a(nu_rest))
    end
end

"""
    build_template_grid(templates::Vector{String}, zgrid::Vector{Float64},
                       igm_model_name::String, fitting_params::Dict; 
                       output_type::Symbol=:photometry, bands::Union{Nothing,Vector{String}}=nothing, 
                       wavelength_grid::Union{Nothing,Vector{Float64}}=nothing)

Unified template grid builder supporting both photometry and spectral output.

Arguments:
- templates: Vector of template file paths
- zgrid: Redshift grid for evaluation
- igm_model_name: Name of IGM model to apply
- fitting_params: Dictionary containing fitting parameters
- output_type: :photometry (for fitting) or :spectral (for template export)
- bands: Required for :photometry mode - filter names
- wavelength_grid: Optional for :spectral mode (auto-determined if nothing)

Returns:
For :photometry mode: (templgrid[nband,ntempl,nz], template_error_grid[nz,nband], nothing)
For :spectral mode: (templgrid[nwav,ntempl,nz], nothing, wavelength_grid)

The template grid shape is (dimension3, ntempl, nz) for both modes, optimized for
column-major access when slicing by redshift index.
"""
function build_template_grid(templates::Vector{String}, zgrid::Vector{Float64},
                           igm_model_name::String, fitting_params::Dict;
                           output_type::Symbol=:photometry, bands::Union{Nothing,Vector{String}}=nothing,
                           wavelength_grid::Union{Nothing,Vector{Float64}}=nothing,
                           add_cgm::Bool=true, cgm_A::Float64=3.5918,
                           cgm_a::Float64=1.8414, cgm_c::Float64=18.001)
    
    # Validate inputs
    if output_type == :photometry && bands === nothing
        throw(LazyError("bands must be provided for photometry mode"))
    end
    
    ntempl = length(templates)
    nz = length(zgrid)
    
    # Determine integration dimension and grid
    if output_type == :photometry
        integration_grid = bands
        nintegration = length(bands)
        println("🧠 Building photometric template grid: ($ntempl, $nz, $nintegration)")
    elseif output_type == :spectral
        # Use provided wavelength grid or determine from templates
        if wavelength_grid === nothing
            wavelength_grid = determine_common_wavelength_grid(templates)
        end
        integration_grid = wavelength_grid
        nintegration = length(wavelength_grid)
        println("🧠 Building spectral template grid: ($ntempl, $nz, $nintegration)")
    else
        throw(LazyError("output_type must be :photometry or :spectral"))
    end
    
    # Initialize template grid: (nintegration, ntempl, nz) for column-major efficiency
    templgrid = zeros(nintegration, ntempl, nz)
    
    # Load IGM model data (same for both modes)
    igm_redshifts, igm_wavelengths, igm_transmission, idx_igm, maxz = load_igm_model(igm_model_name)
    pr = Atomic{Bool}(true)  # Thread-safe print warning flag for high redshift
    
    # Initialize progress bar for template processing
    template_progress = ProgressBar(ntempl, "Processing templates...", 
                                   show_rate=false, show_eta=false, 
                                   prefix_emoji="📊")
    
    # Process each template
    for i in 1:ntempl
        templwav_i, templfnu_i, templz_i = load_template(templates[i])
        idx = searchsortedfirst(templwav_i, 1215.67)
        
        # For spectral mode, determine which wavelength grid to use
        if output_type == :spectral
            # Always use the common wavelength grid for spectral output
            working_wavelength_grid = wavelength_grid
            working_idx = searchsortedfirst(wavelength_grid, 1215.67)
            
            # Pre-interpolate redshift-independent templates to common grid
            if templz_i === nothing && wavelength_grid !== templwav_i
                interp = linear_interpolation(templwav_i, templfnu_i, extrapolation_bc=Flat())
                templfnu_i = interp(wavelength_grid)
            end
        else
            # For photometry mode, use template's native wavelength grid
            working_wavelength_grid = templwav_i
            working_idx = idx
        end
        
        # Process each redshift (parallelized)
        @threads for j in 1:nz
            z = zgrid[j]
            wav_obs = working_wavelength_grid .* (1 + z)
            
            # Apply IGM transmission at this redshift
            if z > maxz
                # Thread-safe test-and-set: only first thread to encounter high-z will print
                if atomic_cas!(pr, true, false)
                    println("⚠️ Redshift grid extends beyond the maximum redshift in the IGM model. Using IGM transmission at z=$maxz")
                end
                iz_up = length(igm_redshifts)
            else
                iz_up = searchsortedfirst(igm_redshifts, z)
            end
            
            transmission = igm_transmission[iz_up,:]
            
            # Interpolate IGM transmission to working wavelength grid
            # Handle edge cases for different wavelength grids
            if working_idx > length(working_wavelength_grid)
                # If working_idx is out of bounds (shouldn't happen), use last valid index
                working_idx = length(working_wavelength_grid)
            end
            
            if working_idx > 1
                interp = linear_interpolation([0.0; igm_wavelengths[1:idx_igm-1]], [0.0; transmission[1:idx_igm-1]], extrapolation_bc=Flat())
                y1 = interp(working_wavelength_grid[1:working_idx-1])
            else
                y1 = Float64[]
            end
            
            if working_idx <= length(working_wavelength_grid)
                interp = linear_interpolation([igm_wavelengths[idx_igm:end]; 1226.0], [transmission[idx_igm:end]; 1.0], extrapolation_bc=Flat())
                y2 = interp(working_wavelength_grid[working_idx:end])
            else
                y2 = Float64[]
            end
            
            transmission_interp = [y1; y2]

            if add_cgm
                apply_cgm_damping_wing!(transmission_interp, working_wavelength_grid, z, cgm_A, cgm_a, cgm_c)
            end

            # Ensure transmission_interp matches working_wavelength_grid length
            if length(transmission_interp) != length(working_wavelength_grid)
                # This shouldn't happen, but handle gracefully
                if length(transmission_interp) < length(working_wavelength_grid)
                    # Pad with ones if too short
                    transmission_interp = [transmission_interp; fill(1.0, length(working_wavelength_grid) - length(transmission_interp))]
                else
                    # Truncate if too long
                    transmission_interp = transmission_interp[1:length(working_wavelength_grid)]
                end
            end
            
            # Apply template redshift dependence if present
            if templz_i === nothing
                # Redshift-independent template
                templfnu_j = templfnu_i .* transmission_interp
            else
                # Redshift-dependent template
                zindex = argmin(abs.(templz_i .- z))
                if output_type == :spectral
                    # For spectral mode, interpolate template to common wavelength grid
                    interp = linear_interpolation(templwav_i, templfnu_i[:,zindex], extrapolation_bc=Flat())
                    templfnu_j = interp(working_wavelength_grid) .* transmission_interp
                else
                    # For photometry mode
                    templfnu_j = templfnu_i[:,zindex] .* transmission_interp
                end
            end
            
            # Normalize template to rest-1micron (numerical stability)
            if output_type == :spectral
                # For spectral mode, normalize based on wavelength_grid
                windex = argmin(abs.(wavelength_grid .- 1e4))
            else
                # For photometry mode, normalize based on template wavelength
                windex = argmin(abs.(templwav_i .- 1e4))
            end
            templfnu_j /= templfnu_j[windex]
            
            # Integration step differs by mode
            if output_type == :photometry
                # Integrate with filter transmission curves
                for k in 1:nintegration
                    band = bands[k]
                    fwav, ftrans = get_filter(band)
                    nu = 1 ./ fwav
                    interp = linear_interpolation(wav_obs, templfnu_j)
                    fnu_interp = interp(fwav)
                    
                    result = trapz(nu, fnu_interp .* ftrans ./ nu) / trapz(nu, ftrans ./ nu)
                    templgrid[k, i, j] = result
                end
                
            elseif output_type == :spectral
                # Sample on wavelength grid (no filter integration)
                templgrid[:, i, j] = templfnu_j
            end
        end
        
        # Update progress
        increment!(template_progress, force=true)
    end
    
    # Finish progress bar
    finish!(template_progress)
    
    # Build template error grid only for photometry mode
    template_error_grid = nothing
    if output_type == :photometry
        if !haskey(fitting_params, "template_error")
            throw(LazyError("parameter `template_error` not found in fitting parameters"))
        end
        if !haskey(fitting_params, "template_error_scale")
            throw(LazyError("parameter `template_error_scale` not found in fitting parameters"))
        end
        
        template_error = fitting_params["template_error"]
        template_error_scale = fitting_params["template_error_scale"]
        
        # Load template error function
        tef = readdlm(templatepath * "template_error/" * template_error * ".dat")
        tef_x = tef[:,1]
        tef_y = tef[:,2]
        tef_xmin = minimum(tef_x[tef_y .> 0])
        tef_xmax = maximum(tef_x[tef_y .> 0])
        tef_clip_min = tef_y[tef_y .> 0][1]
        tef_clip_max = tef_y[tef_y .> 0][end]
        
        # Get pivot wavelengths for bands
        pivot_wavs = get_pivot_wavelengths(bands)
        
        # Pre-compute template error grid
        template_error_grid = precompute_template_errors(zgrid, pivot_wavs, tef_x, tef_y, 
            tef_xmin, tef_xmax, tef_clip_min, tef_clip_max, template_error_scale)
    end
    
    # Return based on mode
    if output_type == :photometry
        return (templgrid, template_error_grid, nothing)
    else  # :spectral
        return (templgrid, nothing, wavelength_grid)
    end
end

"""
    build_restframe_template_grid(templates::Vector{String}, zgrid::Vector{Float64})

Build rest-frame template grid integrated through 4 hardcoded rest-frame filters
(M_UV, M_U, M_V, M_J). No IGM attenuation, no redshifting — just rest-frame F_nu
through rest-frame bandpasses. Normalization matches build_template_grid() (at
rest-frame 10000A).

Returns: restframe_templgrid[4, ntempl, nz] where dim 1 = [M_UV, M_U, M_V, M_J]
"""
function build_restframe_template_grid(templates::Vector{String}, zgrid::Vector{Float64})
    ntempl = length(templates)
    nz = length(zgrid)
    nrf = 4

    # Hardcoded rest-frame filter nicknames
    rf_nicknames = ["rf_uv", "rf_u", "rf_v", "rf_j"]
    rf_filters = [get_filter(nick) for nick in rf_nicknames]

    restframe_templgrid = zeros(nrf, ntempl, nz)

    for i in 1:ntempl
        templwav_i, templfnu_i, templz_i = load_template(templates[i])

        # Normalization index at rest-frame 10000A (same as build_template_grid line 263)
        windex = argmin(abs.(templwav_i .- 1e4))

        @threads for j in 1:nz
            z = zgrid[j]

            # Select flux for this redshift
            if templz_i === nothing
                templfnu_j = copy(templfnu_i)
            else
                zindex = argmin(abs.(templz_i .- z))
                templfnu_j = templfnu_i[:, zindex]
            end

            # Normalize at rest-frame 10000A
            norm_val = templfnu_j[windex]
            if norm_val <= 0 || !isfinite(norm_val)
                continue
            end
            templfnu_j = templfnu_j ./ norm_val

            # Integrate through each rest-frame filter
            for k in 1:nrf
                fwav, ftrans = rf_filters[k]
                nu = 1.0 ./ fwav
                interp = linear_interpolation(templwav_i, templfnu_j, extrapolation_bc=Flat())
                fnu_interp = interp(fwav)

                denom = trapz(nu, ftrans ./ nu)
                if denom != 0
                    restframe_templgrid[k, i, j] = trapz(nu, fnu_interp .* ftrans ./ nu) / denom
                end
            end
        end
    end

    return restframe_templgrid
end