"""
Input/Output operations for Lazy.jl

Contains data loading, template/filter management, HDF5 operations, caching, and FITS I/O.
Merges functionality from main.jl, hdf5_streaming.jl, and cache_utils.jl.
"""

using TOML
using FITSIO
using HDF5
using OrderedCollections
using DelimitedFiles
using Interpolations
using Trapz

# Global paths (from main.jl)
const igmpath = @__DIR__() * "/igm_data/"
const templatepath = @__DIR__() * "/templates/"
const filterpath = @__DIR__() * "/filter_files/"

# Filter caches to avoid repeated disk I/O
const _filter_cache = Dict{String, Tuple{Vector{Float64}, Vector{Float64}}}()
const _filter_directory_cache = Ref{Union{Nothing, Dict}}(nothing)

# =============================================================================
# Data Loading Functions (from main.jl)
# =============================================================================

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
    Load the filters from the filter directory (cached after first call).
    """
    if _filter_directory_cache[] === nothing
        if !isdir(filterpath)
            throw(LazyError("Filter directory $filterpath not found"))
        end
        _filter_directory_cache[] = TOML.parsefile(filterpath * "filter_directory.toml")
    end
    return _filter_directory_cache[]
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

function get_filter(nickname::String)::Tuple{Vector{Float64}, Vector{Float64}}
    """
    Get the filter transmission curve from the nickname (cached after first load).

    Returns (wavelength, transmission) vectors.
    """
    if haskey(_filter_cache, nickname)
        return _filter_cache[nickname]
    end

    filter_directory = load_filters()
    for key in keys(filter_directory)
        for n in filter_directory[key]["nicknames"]
            if n == nickname
                filt = readdlm(filterpath * key)
                result = (filt[:,1], filt[:,2])
                _filter_cache[nickname] = result
                return result
            end
        end
    end
    throw(LazyError("Filter '$nickname' not found"))
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

# =============================================================================
# HDF5 Operations (from hdf5_streaming.jl)
# =============================================================================

"""
    store_param_metadata(parent_group::HDF5.Group, param::Dict)

Store parameter file contents as structured HDF5 groups and datasets.
This avoids issues with TOML.print() returning Nothing.
"""
function store_param_metadata(parent_group::HDF5.Group, param::Dict)
    param_group = create_group(parent_group, "parameters")
    
    for (section_name, section_data) in param
        if isa(section_data, Dict)
            # Create group for each section (io, fitting, runtime, translate)
            section_group = create_group(param_group, section_name)
            
            for (key, value) in section_data
                # Store individual parameters as datasets
                try
                    if isa(value, Dict)
                        # Handle nested dictionaries (like in translate section)
                        nested_group = create_group(section_group, key)
                        for (nested_key, nested_value) in value
                            nested_group[nested_key] = string(nested_value)
                        end
                    elseif isa(value, Vector)
                        # Handle arrays (like template file lists)
                        section_group[key] = string.(value)
                    else
                        # Handle scalar values
                        section_group[key] = string(value)
                    end
                catch e
                    # If storage fails, store as string representation
                    println("Warning: Could not store parameter $section_name.$key as structured data, storing as string")
                    section_group[key] = string(value)
                end
            end
        else
            # Handle top-level non-dict values
            param_group[section_name] = string(section_data)
        end
    end
end

"""
    create_hdf5_work_file(filename::String, param::Dict, nobj::Int, nz::Int, 
                          nband::Int, ntempl::Int, bands::Vector{String}, 
                          zgrid::Vector{Float64}, templates::Vector{String})

Create HDF5 work file with proper structure and metadata for streaming writes.
"""
function create_hdf5_work_file(filename::String, param::Dict, nobj::Int, nz::Int, 
                              nband::Int, ntempl::Int, bands::Vector{String}, 
                              zgrid::Vector{Float64}, templates::Vector{String})
    
    h5open(filename, "w") do file
        # Create metadata group
        g_meta = create_group(file, "metadata")
        
        # Store parameter file as structured HDF5 groups instead of TOML string
        store_param_metadata(g_meta, param)
        
        # Store version safely
        version_str = try
            version
        catch
            "unknown"
        end
        g_meta["version"] = version_str
        
        g_meta["start_time"] = time()
        g_meta["nobj"] = nobj
        g_meta["nz"] = nz
        g_meta["nband"] = nband
        g_meta["ntempl"] = ntempl
        g_meta["bands"] = bands
        g_meta["zgrid"] = zgrid
        g_meta["templates"] = templates
        attrs(g_meta)["last_processed_object"] = 0
        attrs(g_meta)["last_update_time"] = time()
        attrs(g_meta)["complete"] = false
        
        # Create results group with pre-allocated datasets
        g_results = create_group(file, "results")
        create_dataset(g_results, "ID", Int32, (nobj,))
        create_dataset(g_results, "zbest", Float64, (nobj,))
        create_dataset(g_results, "chi2best", Float64, (nobj,))
        create_dataset(g_results, "z_l95", Float64, (nobj,))
        create_dataset(g_results, "z_l68", Float64, (nobj,))
        create_dataset(g_results, "z_med", Float64, (nobj,))
        create_dataset(g_results, "z_u68", Float64, (nobj,))
        create_dataset(g_results, "z_u95", Float64, (nobj,))
        create_dataset(g_results, "coeffs", Float64, (nobj, ntempl))

        # Optional z_spec column
        if get(param["fitting"], "use_zspec", false)
            create_dataset(g_results, "z_spec", Float64, (nobj,))
        end

        # Create photometry group for best-fit model
        g_phot = create_group(file, "photometry")
        for (i, band) in enumerate(bands)
            create_dataset(g_phot, band, Float64, (nobj,))
        end
        
        # Integrated P(z) quantities (always computed)
        z_integers = collect(1:floor(Int, maximum(zgrid)))
        for Z in z_integers
            create_dataset(g_results, "pz_gt$(Z)", Float32, (nobj,))
        end
        create_dataset(g_results, "Sz", Float64, (nobj,))

        # Pz bins (integer-centered unit bins)
        Sz_bc = collect(0.0:1.0:maximum(zgrid))
        for bc in Sz_bc
            create_dataset(g_results, "Pz$(Int(bc))", Float32, (nobj,))
        end
        create_dataset(g_results, "Pz_cen", Float32, (nobj,))
        create_dataset(g_results, "Pz_zgtrzb2", Float32, (nobj,))

        # Optional full P(z) grid with compression
        output_pz = get(param["io"], "output_pz", false)
        if output_pz
            g_pz = create_group(file, "pz")
            # Use chunking and compression for large chi2 grid
            chunk_size = min(1000, nobj)
            create_dataset(g_pz, "chi2grid", Float32, (nz, nobj),
                          chunk=(nz, chunk_size))
            # Store normalized P(z) alongside chi2grid for consistent output
            create_dataset(g_pz, "pz", Float32, (nz, nobj),
                          chunk=(nz, chunk_size))
        end
        
        # Rest-frame absolute magnitudes
        output_restframe_mags = get(param["io"], "output_restframe_mags", false)
        if output_restframe_mags
            for mag_name in ["M_UV", "M_U", "M_V", "M_J"]
                create_dataset(g_results, mag_name, Float64, (nobj,))
            end
        end

        # Forced low-z fit
        output_forced_lowz = get(param["io"], "output_forced_lowz", false)
        if output_forced_lowz
            create_dataset(g_results, "zbest_lowz", Float64, (nobj,))
            create_dataset(g_results, "chi2best_lowz", Float64, (nobj,))
            create_dataset(g_results, "delta_chi2", Float64, (nobj,))
            create_dataset(g_results, "z_l95_lowz", Float64, (nobj,))
            create_dataset(g_results, "z_l68_lowz", Float64, (nobj,))
            create_dataset(g_results, "z_med_lowz", Float64, (nobj,))
            create_dataset(g_results, "z_u68_lowz", Float64, (nobj,))
            create_dataset(g_results, "z_u95_lowz", Float64, (nobj,))
            create_dataset(g_results, "coeffs_lowz", Float64, (nobj, ntempl))
            for band in bands
                create_dataset(g_phot, "$(band)_lowz", Float64, (nobj,))
            end
        end

        # Optional templates group (datasets created when template data is written)
        output_templ = get(param["io"], "output_templates", false)
        if output_templ
            g_templates = create_group(file, "templates")
            # Template datasets will be created dynamically when template data is written
            # This avoids pre-allocating without knowing wavelength grid size
        end
    end
end

"""
    check_resume_file(filename::String)

Check if a work file exists and is resumable. Returns (can_resume, last_object).
"""
function check_resume_file(filename::String)
    if !isfile(filename)
        return (false, 0)
    end
    
    try
        h5open(filename, "r") do file
            if !haskey(file, "metadata")
                return (false, 0)
            end
            
            meta = file["metadata"]
            complete = read_attribute(meta, "complete")
            if complete
                return (false, 0)  # Already complete
            end

            last_obj = read_attribute(meta, "last_processed_object")
            return (true, last_obj)
        end
    catch
        return (false, 0)  # Corrupted file
    end
end

"""
    write_chunk_results(file::HDF5.File, chunk_start, chunk_end, ...)

Write a chunk of results to an open HDF5 file handle.
"""
function write_chunk_results(file::HDF5.File, chunk_start::Int, chunk_end::Int,
                           IDs::Vector{<:Integer}, zbest::Vector{Float64},
                           chi2best::Vector{Float64}, coeffsbest::Matrix{Float64},
                           z_l95::Vector{Float64}, z_l68::Vector{Float64},
                           z_med::Vector{Float64}, z_u68::Vector{Float64},
                           z_u95::Vector{Float64}, photobest::Matrix{Float64},
                           bands::Vector{String}, chi2grid::Union{Nothing,Matrix{Float32}}=nothing,
                           zgrid::Union{Nothing,Vector{Float64}}=nothing;
                           pz_gt::Union{Nothing,Matrix{Float32}}=nothing,
                           Sz::Union{Nothing,Vector{Float64}}=nothing,
                           z_integers::Union{Nothing,Vector{Int}}=nothing,
                           pz_bins::Union{Nothing,Matrix{Float32}}=nothing,
                           pz_bin_centers::Union{Nothing,Vector{Float64}}=nothing,
                           pz_cen::Union{Nothing,Vector{Float32}}=nothing,
                           pz_zgtrzb2::Union{Nothing,Vector{Float32}}=nothing,
                           restframe_mags::Union{Nothing,Matrix{Float64}}=nothing,
                           forced_lowz::Union{Nothing,Dict}=nothing,
                           pz_grid::Union{Nothing,Matrix{Float32}}=nothing)

    # Write results (convert IDs to Int32 to match dataset type)
    file["results/ID"][chunk_start:chunk_end] = Int32.(IDs)
    file["results/zbest"][chunk_start:chunk_end] = zbest
    file["results/chi2best"][chunk_start:chunk_end] = chi2best
    file["results/z_l95"][chunk_start:chunk_end] = z_l95
    file["results/z_l68"][chunk_start:chunk_end] = z_l68
    file["results/z_med"][chunk_start:chunk_end] = z_med
    file["results/z_u68"][chunk_start:chunk_end] = z_u68
    file["results/z_u95"][chunk_start:chunk_end] = z_u95
    file["results/coeffs"][chunk_start:chunk_end, :] = coeffsbest

    # Write photometry
    for (i, band) in enumerate(bands)
        file["photometry/$band"][chunk_start:chunk_end] = photobest[:, i]
    end

    # Write P(z) if provided
    if chi2grid !== nothing && haskey(file, "pz/chi2grid")
        file["pz/chi2grid"][:, chunk_start:chunk_end] = chi2grid

        # Write pre-computed normalized P(z) if provided, otherwise compute it
        if haskey(file, "pz/pz")
            if pz_grid !== nothing
                file["pz/pz"][:, chunk_start:chunk_end] = pz_grid
            elseif zgrid !== nothing
                pz_chunk = exp.(-0.5f0 .* chi2grid)
                for col in 1:size(pz_chunk, 2)
                    norm = trapz(zgrid, @view pz_chunk[:, col])
                    if norm > 0
                        pz_chunk[:, col] ./= norm
                    end
                end
                file["pz/pz"][:, chunk_start:chunk_end] = pz_chunk
            end
        end
    end

    # Write integrated P(z) quantities
    if pz_gt !== nothing && z_integers !== nothing
        for (zi, Z) in enumerate(z_integers)
            if haskey(file["results"], "pz_gt$(Z)")
                file["results/pz_gt$(Z)"][chunk_start:chunk_end] = pz_gt[:, zi]
            end
        end
    end
    if Sz !== nothing && haskey(file["results"], "Sz")
        file["results/Sz"][chunk_start:chunk_end] = Sz
    end

    # Write Pz bin probabilities
    if pz_bins !== nothing && pz_bin_centers !== nothing
        for (i, bc) in enumerate(pz_bin_centers)
            key = "Pz$(Int(bc))"
            if haskey(file["results"], key)
                file["results/$key"][chunk_start:chunk_end] = pz_bins[:, i]
            end
        end
    end
    if pz_cen !== nothing && haskey(file["results"], "Pz_cen")
        file["results/Pz_cen"][chunk_start:chunk_end] = pz_cen
    end
    if pz_zgtrzb2 !== nothing && haskey(file["results"], "Pz_zgtrzb2")
        file["results/Pz_zgtrzb2"][chunk_start:chunk_end] = pz_zgtrzb2
    end

    # Write rest-frame absolute magnitudes
    if restframe_mags !== nothing
        for (k, mag_name) in enumerate(["M_UV", "M_U", "M_V", "M_J"])
            if haskey(file["results"], mag_name)
                file["results/$mag_name"][chunk_start:chunk_end] = restframe_mags[:, k]
            end
        end
    end

    # Write forced low-z results
    if forced_lowz !== nothing
        if haskey(file["results"], "zbest_lowz")
            file["results/zbest_lowz"][chunk_start:chunk_end] = forced_lowz["zbest"]
            file["results/chi2best_lowz"][chunk_start:chunk_end] = forced_lowz["chi2best"]
            file["results/delta_chi2"][chunk_start:chunk_end] = forced_lowz["delta_chi2"]
            file["results/z_l95_lowz"][chunk_start:chunk_end] = forced_lowz["z_l95"]
            file["results/z_l68_lowz"][chunk_start:chunk_end] = forced_lowz["z_l68"]
            file["results/z_med_lowz"][chunk_start:chunk_end] = forced_lowz["z_med"]
            file["results/z_u68_lowz"][chunk_start:chunk_end] = forced_lowz["z_u68"]
            file["results/z_u95_lowz"][chunk_start:chunk_end] = forced_lowz["z_u95"]
            file["results/coeffs_lowz"][chunk_start:chunk_end, :] = forced_lowz["coeffs"]
            for (i, band) in enumerate(bands)
                if haskey(file["photometry"], "$(band)_lowz")
                    file["photometry/$(band)_lowz"][chunk_start:chunk_end] = forced_lowz["photobest"][:, i]
                end
            end
        end
    end

    # Update metadata
    attrs(file["metadata"])["last_processed_object"] = chunk_end
    attrs(file["metadata"])["last_update_time"] = time()
end

"""
    write_chunk_results(filename::String, chunk_start, chunk_end, ...)

Convenience wrapper that opens the HDF5 file, writes, and flushes.
"""
function write_chunk_results(filename::String, chunk_start::Int, chunk_end::Int,
                           IDs::Vector{<:Integer}, zbest::Vector{Float64},
                           chi2best::Vector{Float64}, coeffsbest::Matrix{Float64},
                           z_l95::Vector{Float64}, z_l68::Vector{Float64},
                           z_med::Vector{Float64}, z_u68::Vector{Float64},
                           z_u95::Vector{Float64}, photobest::Matrix{Float64},
                           bands::Vector{String}, chi2grid::Union{Nothing,Matrix{Float32}}=nothing,
                           zgrid::Union{Nothing,Vector{Float64}}=nothing;
                           pz_gt::Union{Nothing,Matrix{Float32}}=nothing,
                           Sz::Union{Nothing,Vector{Float64}}=nothing,
                           z_integers::Union{Nothing,Vector{Int}}=nothing,
                           pz_bins::Union{Nothing,Matrix{Float32}}=nothing,
                           pz_bin_centers::Union{Nothing,Vector{Float64}}=nothing,
                           pz_cen::Union{Nothing,Vector{Float32}}=nothing,
                           pz_zgtrzb2::Union{Nothing,Vector{Float32}}=nothing,
                           restframe_mags::Union{Nothing,Matrix{Float64}}=nothing,
                           forced_lowz::Union{Nothing,Dict}=nothing,
                           pz_grid::Union{Nothing,Matrix{Float32}}=nothing)
    try
        h5open(filename, "r+") do file
            write_chunk_results(file, chunk_start, chunk_end,
                              IDs, zbest, chi2best, coeffsbest,
                              z_l95, z_l68, z_med, z_u68, z_u95,
                              photobest, bands, chi2grid, zgrid;
                              pz_gt=pz_gt, Sz=Sz, z_integers=z_integers,
                              pz_bins=pz_bins, pz_bin_centers=pz_bin_centers,
                              pz_cen=pz_cen, pz_zgtrzb2=pz_zgtrzb2,
                              restframe_mags=restframe_mags, forced_lowz=forced_lowz,
                              pz_grid=pz_grid)
            flush(file)
        end
    catch e
        println("Error writing chunk results to $filename:")
        println("   Chunk: objects $chunk_start to $chunk_end")
        println("   Array sizes: IDs=$(length(IDs)), zbest=$(length(zbest))")
        error_msg = if hasfield(typeof(e), :msg)
            e.msg
        else
            string(e)
        end
        println("   Error: $(typeof(e)) - $error_msg")
        rethrow(e)
    end
end

"""
    write_template_data(work_file::String, templgrid_spectral, wavelength_grid, zgrid, templates)

Write template spectral data to HDF5 work file.
templgrid_spectral has shape (nwav, ntempl, nz). Per template, slice [:, i, :] gives (nwav, nz).
"""
function write_template_data(work_file::String, templgrid_spectral, wavelength_grid, zgrid, templates)
    h5open(work_file, "r+") do file
        g_templates = file["templates"]

        # Store wavelength grid (common to all templates)
        g_templates["wavelength"] = wavelength_grid

        # Store each template's spectral data
        # templgrid_spectral[:, i, :] has shape (nwav, nz) — already correct orientation
        for (i, template) in enumerate(templates)
            template_name = basename(template)
            template_name = splitext(template_name)[1]  # Remove extension regardless of type

            g_templates[template_name] = templgrid_spectral[:, i, :]
        end
    end
end

"""
    finalize_hdf5_work_file(filename::String)

Mark HDF5 work file as complete.
"""
function finalize_hdf5_work_file(filename::String)
    h5open(filename, "r+") do file
        meta_group = file["metadata"]
        attrs(meta_group)["complete"] = true
        attrs(meta_group)["end_time"] = time()
    end
end

"""
    convert_hdf5_to_fits_python(hdf5_file::String, fits_file::String)

Convert completed HDF5 work file to FITS format using PyCall/astropy.
Preserved as fallback; use convert_hdf5_to_fits() for the default CFITSIO-based writer.
"""
function convert_hdf5_to_fits_python(hdf5_file::String, fits_file::String)
    h5open(hdf5_file, "r") do h5f
        # Read metadata
        meta = h5f["metadata"]
        bands = read(meta, "bands")
        nobj = read(meta, "nobj")

        # Read results
        IDs = read(h5f["results/ID"])
        zbest = read(h5f["results/zbest"])
        chi2best = read(h5f["results/chi2best"])
        z_l95 = read(h5f["results/z_l95"])
        z_l68 = read(h5f["results/z_l68"])
        z_med = read(h5f["results/z_med"])
        z_u68 = read(h5f["results/z_u68"])
        z_u95 = read(h5f["results/z_u95"])
        coeffsbest = read(h5f["results/coeffs"])

        # Prepare data for FITS
        data = OrderedDict{String, Dict{String, Any}}()
        data["ID"] = Dict("format" => "K", "data" => IDs)
        data["z_best"] = Dict("format" => "E", "data" => zbest)
        data["chi2"] = Dict("format" => "E", "data" => chi2best)
        # Include z_spec if available
        if haskey(h5f["results"], "z_spec")
            z_spec_data = read(h5f["results/z_spec"])
            data["z_spec"] = Dict("format" => "E", "data" => z_spec_data)
        end
        data["z_l95"] = Dict("format" => "E", "data" => z_l95)
        data["z_l68"] = Dict("format" => "E", "data" => z_l68)
        data["z_med"] = Dict("format" => "E", "data" => z_med)
        data["z_u68"] = Dict("format" => "E", "data" => z_u68)
        data["z_u95"] = Dict("format" => "E", "data" => z_u95)

        # Integrated P(z) columns
        zgrid_py = read(meta, "zgrid")
        z_integers_py = collect(1:floor(Int, maximum(zgrid_py)))
        has_pz_gt_py = !isempty(z_integers_py) && haskey(h5f["results"], "pz_gt$(z_integers_py[1])")
        if has_pz_gt_py
            for Z in z_integers_py
                data["Pz_gt$(Z)"] = Dict("format" => "E", "data" => read(h5f["results/pz_gt$(Z)"]))
            end
        end
        has_Sz_py = haskey(h5f["results"], "Sz")
        if has_Sz_py
            data["Sz"] = Dict("format" => "E", "data" => read(h5f["results/Sz"]))
        end

        # Pz bin columns
        Sz_bc_py = collect(0:floor(Int, maximum(zgrid_py)))
        has_pz_bins_py = !isempty(Sz_bc_py) && haskey(h5f["results"], "Pz$(Sz_bc_py[1])")
        if has_pz_bins_py
            for bc in Sz_bc_py
                data["Pz$(bc)"] = Dict("format" => "E", "data" => read(h5f["results/Pz$(bc)"]))
            end
        end
        if haskey(h5f["results"], "Pz_cen")
            data["Pz_cen"] = Dict("format" => "E", "data" => read(h5f["results/Pz_cen"]))
        end
        if haskey(h5f["results"], "Pz_zgtrzb2")
            data["Pz_zgtrzb2"] = Dict("format" => "E", "data" => read(h5f["results/Pz_zgtrzb2"]))
        end

        # Rest-frame magnitudes
        has_restframe_mags_py = haskey(h5f["results"], "M_UV")
        if has_restframe_mags_py
            for mag_name in ["M_UV", "M_U", "M_V", "M_J"]
                data[mag_name] = Dict("format" => "E", "unit" => "mag", "data" => read(h5f["results/$mag_name"]))
            end
        end

        # Forced low-z columns
        has_forced_lowz_py = haskey(h5f["results"], "zbest_lowz")
        if has_forced_lowz_py
            data["z_best_lowz"] = Dict("format" => "E", "data" => read(h5f["results/zbest_lowz"]))
            data["chi2_lowz"] = Dict("format" => "E", "data" => read(h5f["results/chi2best_lowz"]))
            data["delta_chi2"] = Dict("format" => "E", "data" => read(h5f["results/delta_chi2"]))
            data["z_l95_lowz"] = Dict("format" => "E", "data" => read(h5f["results/z_l95_lowz"]))
            data["z_l68_lowz"] = Dict("format" => "E", "data" => read(h5f["results/z_l68_lowz"]))
            data["z_med_lowz"] = Dict("format" => "E", "data" => read(h5f["results/z_med_lowz"]))
            data["z_u68_lowz"] = Dict("format" => "E", "data" => read(h5f["results/z_u68_lowz"]))
            data["z_u95_lowz"] = Dict("format" => "E", "data" => read(h5f["results/z_u95_lowz"]))
        end

        # Add photometry
        for (i, band) in enumerate(bands)
            photobest_band = read(h5f["photometry/$band"])
            data[band] = Dict("format" => "E", "unit" => "fnu", "data" => photobest_band)
        end

        ntempl = size(coeffsbest, 2)
        data["coeffs"] = Dict("format" => "$(ntempl)E", "data" => coeffsbest)

        if has_forced_lowz_py
            for (i, band) in enumerate(bands)
                data["$(band)_lowz"] = Dict("format" => "E", "unit" => "fnu", "data" => read(h5f["photometry/$(band)_lowz"]))
            end
            coeffs_lowz_py = read(h5f["results/coeffs_lowz"])
            data["coeffs_lowz"] = Dict("format" => "$(ntempl)E", "data" => coeffs_lowz_py)
        end

        # Write SUMMARY extension
        write_data(fits_file, data, "SUMMARY")

        # Write P(z) if it exists
        if haskey(h5f, "pz/chi2grid")
            zgrid = read(meta, "zgrid")
            chi2grid = read(h5f["pz/chi2grid"])  # (nz, nobj)

            # Convert chi2 to P(z) and normalize each column (object)
            pz = exp.(-0.5 * chi2grid)
            for col in 1:size(pz, 2)
                norm = trapz(zgrid, @view pz[:, col])
                if norm > 0
                    pz[:, col] ./= norm
                end
            end

            # Format for FITS: (nobj+1, nz) with zgrid as first row
            temp_pz = vcat(transpose(zgrid), transpose(pz))
            temp_IDs = vcat([-1], IDs)

            nz = length(zgrid)
            data_pz = OrderedDict{String, Dict{String, Any}}()
            data_pz["ID"] = Dict("format" => "K", "data" => temp_IDs)
            data_pz["Pz"] = Dict("format" => "$(nz)E", "data" => temp_pz)

            write_data(fits_file, data_pz, "PZ")
        end

        # Write TEMPL extension if it exists
        if haskey(h5f, "templates")
            zgrid = read(meta, "zgrid")
            templates = read(meta, "templates")
            wavelength_grid = read(h5f["templates/wavelength"])

            nwav = length(wavelength_grid)

            # Prepare data for FITS (matching original format exactly)
            data_templ = OrderedDict{String, Dict{String, Any}}()
            temp_zgrid = vcat([-1], zgrid)
            data_templ["z"] = Dict("format" => "E", "data" => temp_zgrid)

            for template in templates
                template_name = basename(template)
                template_name = splitext(template_name)[1]

                # Read template data (stored as (nwav, nz) format)
                template_data = read(h5f["templates/$template_name"])

                # Format for FITS output (matching original exactly)
                # Need to transpose template_data back to (nz, nwav) for FITS format
                # FITS format: (nz+1, nwav) where first row is wavelength, then nz rows of template flux
                template_data_transposed = transpose(template_data)  # Convert (nwav, nz) -> (nz, nwav)
                temp_template_data = vcat(transpose(wavelength_grid), template_data_transposed)
                data_templ[template_name] = Dict("format" => "$(nwav)E", "data" => temp_template_data)
            end

            write_data(fits_file, data_templ, "TEMPL")
        end
    end
end

"""
    convert_hdf5_to_fits(hdf5_file::String, fits_file::String; chunk_size::Int=100_000)

Convert completed HDF5 work file to FITS format using CFITSIO.jl with chunked writes.
Reads data from HDF5 in chunks to avoid loading everything into memory at once.
"""
function convert_hdf5_to_fits(hdf5_file::String, fits_file::String; chunk_size::Int=100_000)
    h5open(hdf5_file, "r") do h5f
        meta = h5f["metadata"]
        nobj = read(meta, "nobj")
        bands = read(meta, "bands")
        ntempl = read(meta, "ntempl")
        has_zspec = haskey(h5f["results"], "z_spec")

        # Remove existing file if present
        isfile(fits_file) && rm(fits_file)

        ff = CFITSIO.fits_clobber_file(fits_file)
        try
            # ── SUMMARY extension ──
            # Build column definitions: (name, tform, unit)
            summary_coldefs = NTuple{3,String}[
                ("ID", "1K", ""),
                ("z_best", "1E", ""),
                ("chi2", "1E", ""),
            ]
            if has_zspec
                push!(summary_coldefs, ("z_spec", "1E", ""))
            end
            append!(summary_coldefs, [
                ("z_l95", "1E", ""),
                ("z_l68", "1E", ""),
                ("z_med", "1E", ""),
                ("z_u68", "1E", ""),
                ("z_u95", "1E", ""),
            ])

            # Integrated P(z) columns
            zgrid_full = read(meta, "zgrid")
            z_integers = collect(1:floor(Int, maximum(zgrid_full)))
            has_pz_gt = !isempty(z_integers) && haskey(h5f["results"], "pz_gt$(z_integers[1])")
            if has_pz_gt
                for Z in z_integers
                    push!(summary_coldefs, ("Pz_gt$(Z)", "1E", ""))
                end
            end
            has_Sz = haskey(h5f["results"], "Sz")
            if has_Sz
                push!(summary_coldefs, ("Sz", "1E", ""))
            end

            # Pz bin columns (integer-centered unit bins)
            Sz_bc_fits = collect(0:floor(Int, maximum(zgrid_full)))
            has_pz_bins = !isempty(Sz_bc_fits) && haskey(h5f["results"], "Pz$(Sz_bc_fits[1])")
            if has_pz_bins
                for bc in Sz_bc_fits
                    push!(summary_coldefs, ("Pz$(bc)", "1E", ""))
                end
            end
            has_pz_cen = haskey(h5f["results"], "Pz_cen")
            if has_pz_cen
                push!(summary_coldefs, ("Pz_cen", "1E", ""))
            end
            has_pz_zgtrzb2 = haskey(h5f["results"], "Pz_zgtrzb2")
            if has_pz_zgtrzb2
                push!(summary_coldefs, ("Pz_zgtrzb2", "1E", ""))
            end

            # Rest-frame magnitude columns
            has_restframe_mags = haskey(h5f["results"], "M_UV")
            if has_restframe_mags
                for mag_name in ["M_UV", "M_U", "M_V", "M_J"]
                    push!(summary_coldefs, (mag_name, "1E", "mag"))
                end
            end

            # Forced low-z columns
            has_forced_lowz = haskey(h5f["results"], "zbest_lowz")
            if has_forced_lowz
                append!(summary_coldefs, [
                    ("z_best_lowz", "1E", ""),
                    ("chi2_lowz", "1E", ""),
                    ("delta_chi2", "1E", ""),
                    ("z_l95_lowz", "1E", ""),
                    ("z_l68_lowz", "1E", ""),
                    ("z_med_lowz", "1E", ""),
                    ("z_u68_lowz", "1E", ""),
                    ("z_u95_lowz", "1E", ""),
                ])
            end

            for band in bands
                push!(summary_coldefs, (band, "1E", "fnu"))
            end
            push!(summary_coldefs, ("coeffs", "$(ntempl)E", ""))

            if has_forced_lowz
                for band in bands
                    push!(summary_coldefs, ("$(band)_lowz", "1E", "fnu"))
                end
                push!(summary_coldefs, ("coeffs_lowz", "$(ntempl)E", ""))
            end

            CFITSIO.fits_create_binary_tbl(ff, nobj, summary_coldefs, "SUMMARY")

            # Write SUMMARY in chunks
            for cs in 1:chunk_size:nobj
                ce = min(cs + chunk_size - 1, nobj)

                colnum = 1
                CFITSIO.fits_write_col(ff, colnum, cs, 1, Int64.(h5f["results/ID"][cs:ce])); colnum += 1
                CFITSIO.fits_write_col(ff, colnum, cs, 1, Float32.(h5f["results/zbest"][cs:ce])); colnum += 1
                CFITSIO.fits_write_col(ff, colnum, cs, 1, Float32.(h5f["results/chi2best"][cs:ce])); colnum += 1
                if has_zspec
                    CFITSIO.fits_write_col(ff, colnum, cs, 1, Float32.(h5f["results/z_spec"][cs:ce])); colnum += 1
                end
                CFITSIO.fits_write_col(ff, colnum, cs, 1, Float32.(h5f["results/z_l95"][cs:ce])); colnum += 1
                CFITSIO.fits_write_col(ff, colnum, cs, 1, Float32.(h5f["results/z_l68"][cs:ce])); colnum += 1
                CFITSIO.fits_write_col(ff, colnum, cs, 1, Float32.(h5f["results/z_med"][cs:ce])); colnum += 1
                CFITSIO.fits_write_col(ff, colnum, cs, 1, Float32.(h5f["results/z_u68"][cs:ce])); colnum += 1
                CFITSIO.fits_write_col(ff, colnum, cs, 1, Float32.(h5f["results/z_u95"][cs:ce])); colnum += 1
                if has_pz_gt
                    for Z in z_integers
                        CFITSIO.fits_write_col(ff, colnum, cs, 1, Float32.(h5f["results/pz_gt$(Z)"][cs:ce])); colnum += 1
                    end
                end
                if has_Sz
                    CFITSIO.fits_write_col(ff, colnum, cs, 1, Float32.(h5f["results/Sz"][cs:ce])); colnum += 1
                end
                if has_pz_bins
                    for bc in Sz_bc_fits
                        CFITSIO.fits_write_col(ff, colnum, cs, 1, Float32.(h5f["results/Pz$(bc)"][cs:ce])); colnum += 1
                    end
                end
                if has_pz_cen
                    CFITSIO.fits_write_col(ff, colnum, cs, 1, Float32.(h5f["results/Pz_cen"][cs:ce])); colnum += 1
                end
                if has_pz_zgtrzb2
                    CFITSIO.fits_write_col(ff, colnum, cs, 1, Float32.(h5f["results/Pz_zgtrzb2"][cs:ce])); colnum += 1
                end
                if has_restframe_mags
                    for mag_name in ["M_UV", "M_U", "M_V", "M_J"]
                        CFITSIO.fits_write_col(ff, colnum, cs, 1, Float32.(h5f["results/$mag_name"][cs:ce])); colnum += 1
                    end
                end
                if has_forced_lowz
                    CFITSIO.fits_write_col(ff, colnum, cs, 1, Float32.(h5f["results/zbest_lowz"][cs:ce])); colnum += 1
                    CFITSIO.fits_write_col(ff, colnum, cs, 1, Float32.(h5f["results/chi2best_lowz"][cs:ce])); colnum += 1
                    CFITSIO.fits_write_col(ff, colnum, cs, 1, Float32.(h5f["results/delta_chi2"][cs:ce])); colnum += 1
                    CFITSIO.fits_write_col(ff, colnum, cs, 1, Float32.(h5f["results/z_l95_lowz"][cs:ce])); colnum += 1
                    CFITSIO.fits_write_col(ff, colnum, cs, 1, Float32.(h5f["results/z_l68_lowz"][cs:ce])); colnum += 1
                    CFITSIO.fits_write_col(ff, colnum, cs, 1, Float32.(h5f["results/z_med_lowz"][cs:ce])); colnum += 1
                    CFITSIO.fits_write_col(ff, colnum, cs, 1, Float32.(h5f["results/z_u68_lowz"][cs:ce])); colnum += 1
                    CFITSIO.fits_write_col(ff, colnum, cs, 1, Float32.(h5f["results/z_u95_lowz"][cs:ce])); colnum += 1
                end
                for band in bands
                    CFITSIO.fits_write_col(ff, colnum, cs, 1, Float32.(h5f["photometry/$band"][cs:ce])); colnum += 1
                end
                # Coefficients: (n, ntempl) matrix — flatten row-major for FITS
                coeffs_chunk = Float32.(h5f["results/coeffs"][cs:ce, :])
                CFITSIO.fits_write_col(ff, colnum, cs, 1, vec(permutedims(coeffs_chunk))); colnum += 1
                if has_forced_lowz
                    for band in bands
                        CFITSIO.fits_write_col(ff, colnum, cs, 1, Float32.(h5f["photometry/$(band)_lowz"][cs:ce])); colnum += 1
                    end
                    coeffs_lowz_chunk = Float32.(h5f["results/coeffs_lowz"][cs:ce, :])
                    CFITSIO.fits_write_col(ff, colnum, cs, 1, vec(permutedims(coeffs_lowz_chunk))); colnum += 1
                end
            end

            # ── PZ extension ──
            if haskey(h5f, "pz/pz")
                zgrid = read(meta, "zgrid")
                nz = length(zgrid)

                pz_coldefs = NTuple{3,String}[
                    ("ID", "1K", ""),
                    ("Pz", "$(nz)E", ""),
                ]
                CFITSIO.fits_create_binary_tbl(ff, nobj + 1, pz_coldefs, "PZ")

                # Row 1: sentinel row with ID=-1 and Pz=zgrid
                CFITSIO.fits_write_col(ff, 1, 1, 1, Int64[-1])
                CFITSIO.fits_write_col(ff, 2, 1, 1, Float32.(zgrid))

                # Remaining rows: chunked from pre-computed pz dataset
                for cs in 1:chunk_size:nobj
                    ce = min(cs + chunk_size - 1, nobj)
                    ids_chunk = Int64.(h5f["results/ID"][cs:ce])
                    pz_chunk = Float32.(h5f["pz/pz"][:, cs:ce])  # (nz, chunk_nobj)

                    # Write to rows cs+1:ce+1 (offset by 1 for sentinel row)
                    CFITSIO.fits_write_col(ff, 1, cs + 1, 1, ids_chunk)
                    CFITSIO.fits_write_col(ff, 2, cs + 1, 1, vec(pz_chunk))
                end
            end

            # ── TEMPL extension ──
            if haskey(h5f, "templates")
                zgrid = read(meta, "zgrid")
                nz = length(zgrid)
                template_names_raw = read(meta, "templates")
                wavelength_grid = read(h5f["templates/wavelength"])
                nwav = length(wavelength_grid)

                templ_coldefs = NTuple{3,String}[("z", "1E", "")]
                template_names = String[]
                for template in template_names_raw
                    tname = splitext(basename(template))[1]
                    push!(template_names, tname)
                    push!(templ_coldefs, (tname, "$(nwav)E", ""))
                end

                nrows_templ = nz + 1  # 1 sentinel row + nz redshift rows
                CFITSIO.fits_create_binary_tbl(ff, nrows_templ, templ_coldefs, "TEMPL")

                # Row 1: sentinel with z=-1, template columns = wavelength grid
                CFITSIO.fits_write_col(ff, 1, 1, 1, Float32[-1])
                for (ti, tname) in enumerate(template_names)
                    CFITSIO.fits_write_col(ff, 1 + ti, 1, 1, Float32.(wavelength_grid))
                end

                # Rows 2..nz+1: redshift grid values and template fluxes
                CFITSIO.fits_write_col(ff, 1, 2, 1, Float32.(zgrid))
                for (ti, tname) in enumerate(template_names)
                    # Template data stored as (nwav, nz) in HDF5; need (nz, nwav) for FITS rows
                    template_data = read(h5f["templates/$tname"])  # (nwav, nz)
                    template_transposed = Float32.(permutedims(template_data))  # (nz, nwav)
                    CFITSIO.fits_write_col(ff, 1 + ti, 2, 1, vec(template_transposed))
                end
            end

        finally
            CFITSIO.fits_close_file(ff)
        end
    end
end

function get_optimal_chunk_size(nobj::Int, nz::Int)
    """
    Calculate optimal chunk size based on default memory target (0.5 GB).
    Maintained for backwards compatibility.
    """
    return get_chunk_size_for_memory_target(nobj, nz, 0.5)
end

"""
    get_chunk_size_for_memory_target(nobj::Int, nz::Int, target_memory_gb::Float64)

Calculate chunk size to achieve target memory usage per chunk.
"""
function get_chunk_size_for_memory_target(nobj::Int, nz::Int, target_memory_gb::Float64)
    bytes_per_object = nz * 4  # Float32 chi2 values
    target_bytes = target_memory_gb * 1024^3
    objects_per_chunk = floor(Int, target_bytes / bytes_per_object)
    
    # Clamp to reasonable range (min 1000 objects, max total objects)
    return clamp(objects_per_chunk, 1_000, nobj)
end

"""
    prompt_complete_run(work_file::String, last_obj::Int, nobj::Int)

Prompt user when a run is already complete.
Returns: action choice ("complete", "overwrite", "keep")
"""
function prompt_complete_run(work_file::String, last_obj::Int, nobj::Int)::String
    # Get file modification time
    mtime = stat(work_file).mtime
    time_ago = time() - mtime
    time_str = if time_ago < 3600
        "$(round(Int, time_ago / 60)) minutes ago"
    elseif time_ago < 86400
        "$(round(time_ago / 3600, digits=1)) hours ago"
    else
        "$(round(Int, time_ago / 86400)) days ago"
    end
    
    println("🧠 Found complete run (100.0% complete, $(format_number(last_obj))/$(format_number(nobj)) objects)")
    println("   Last updated: $time_str")
    println("   Options:")
    println("   [C]omplete - Mark as complete and create final output")
    println("   [O]verwrite - Delete work file and start fresh")
    println("   [K]eep - Keep work file and exit")
    print("   Choice [C/o/k]: ")
    
    response = lowercase(strip(readline()))
    
    if response == "o"
        return "overwrite"
    elseif response == "k"
        return "keep"
    else  # Default to complete
        return "complete"
    end
end

"""
    prompt_resume(work_file::String, last_obj::Int, nobj::Int)

Prompt user whether to resume from checkpoint.
Returns: (action, should_resume) where action is "resume", "overwrite", "complete", or "keep"
"""
function prompt_resume(work_file::String, last_obj::Int, nobj::Int)::Tuple{String, Bool}
    # Check if the run is actually complete
    if last_obj >= nobj
        action = prompt_complete_run(work_file, last_obj, nobj)
        return (action, false)  # No resuming needed for complete runs
    end
    
    percent_complete = round(100 * last_obj / nobj, digits=1)
    
    # Get file modification time
    mtime = stat(work_file).mtime
    time_ago = time() - mtime
    time_str = if time_ago < 3600
        "$(round(Int, time_ago / 60)) minutes ago"
    elseif time_ago < 86400
        "$(round(time_ago / 3600, digits=1)) hours ago"
    else
        "$(round(Int, time_ago / 86400)) days ago"
    end
    
    println("🧠 Found incomplete run ($percent_complete% complete, $(format_number(last_obj))/$(format_number(nobj)) objects)")
    println("   Last updated: $time_str")
    print("   Resume from object $(last_obj + 1)? [Y/n]: ")
    
    response = readline()
    if response == "" || lowercase(response) == "y"
        return ("resume", true)
    else
        return ("overwrite", false)
    end
end

# =============================================================================
# Template Caching Functions (from cache_utils.jl)
# =============================================================================

const CACHE_DIR = joinpath(homedir(), ".lazy", "cache")
const MAX_CACHE_SIZE_GB = 10.0  # Maximum cache size in GB

"""
Ensure the cache directory exists
"""
function ensure_cache_dir()
    if !isdir(CACHE_DIR)
        mkpath(CACHE_DIR)
        println("📁 Created cache directory: $CACHE_DIR")
    end
end

"""
Generate a unique cache key based on template grid parameters
"""
function generate_cache_key(templates::Vector{String}, zgrid::Vector{Float64},
                           bands::Vector{String}, igm_model::String,
                           template_error::String, template_error_scale::Float64;
                           add_cgm::Bool=true, cgm_A::Float64=3.5918,
                           cgm_a::Float64=1.8414, cgm_c::Float64=18.001)::String

    # Create a hash of all relevant parameters
    template_names = sort([basename(t) for t in templates])
    sorted_bands = sort(bands)

    # Include a code version fingerprint so cache is invalidated when
    # template_grid.jl or io.jl change (prevents loading stale grids)
    code_files = [joinpath(@__DIR__, "template_grid.jl"), joinpath(@__DIR__, "io.jl")]
    code_fingerprint = string(sum(stat(f).mtime for f in code_files if isfile(f)))

    # Create a string representation of all parameters
    params_str = join([
        join(template_names, ","),
        "z$(minimum(zgrid))_$(maximum(zgrid))_$(length(zgrid))",
        join(sorted_bands, ","),
        igm_model,
        template_error,
        string(template_error_scale),
        string(add_cgm),
        string(cgm_A), string(cgm_a), string(cgm_c),
        code_fingerprint
    ], "_")
    
    # Generate hash and take first 8 characters for readability
    hash_val = hash(params_str)
    return string(hash_val, base=16)[1:8]
end

"""
Get the full path for a cache file
"""
function get_cache_filename(cache_key::String)::String
    return joinpath(CACHE_DIR, "templgrid_$(cache_key).h5")
end

"""
Check if a cache file exists for the given key
"""
function cache_exists(cache_key::String)::Bool
    filename = get_cache_filename(cache_key)
    return isfile(filename)
end

"""
Save template grid and error grid to cache
"""
function save_template_cache(templgrid::Array{Float32,3}, template_error_grid::Matrix{Float64},
                            metadata::Dict, cache_key::String)
    try
        ensure_cache_dir()
        filename = get_cache_filename(cache_key)
        
        h5open(filename, "w") do file
            # Store the template grid
            file["templgrid"] = templgrid
            file["template_error_grid"] = template_error_grid
            
            # Store metadata
            meta_group = create_group(file, "metadata")
            for (key, value) in metadata
                try
                    if isa(value, Vector{String})
                        meta_group[key] = value
                    elseif isa(value, Vector)
                        meta_group[key] = collect(value)  # Convert to regular array
                    else
                        meta_group[key] = value
                    end
                catch e
                    # For complex types, store as string
                    meta_group[key] = string(value)
                end
            end
            
            # Store cache creation info
            meta_group["cache_created"] = time()
            meta_group["cache_key"] = cache_key
        end
        
        # Clean up old cache files if needed
        cleanup_cache()
        
    catch e
        @warn "Failed to save template cache: $e"
        # Don't let cache failures stop the main computation
    end
end

"""
Load template grid from cache if it exists
"""
function load_template_cache(cache_key::String)::Union{Tuple{Array{Float32,3}, Matrix{Float64}}, Nothing}
    try
        filename = get_cache_filename(cache_key)
        
        if !isfile(filename)
            return nothing
        end
        
        h5open(filename, "r") do file
            templgrid = Float32.(read(file, "templgrid"))
            template_error_grid = read(file, "template_error_grid")
            # Validate shape: expected (nband, ntempl, nz) where dim1 == nband
            nband_cached = size(template_error_grid, 2)
            if size(templgrid, 1) != nband_cached
                @warn "Cached template grid has incompatible shape $(size(templgrid)), rebuilding."
                return nothing
            end
            return (templgrid, template_error_grid)
        end
        
    catch e
        @warn "Failed to load template cache $cache_key: $e"
        return nothing
    end
end

"""
Get information about cache usage
"""
function get_cache_info()::Dict
    if !isdir(CACHE_DIR)
        return Dict("exists" => false)
    end
    
    cache_files = filter(f -> endswith(f, ".h5"), readdir(CACHE_DIR))
    total_size = 0
    file_info = []
    
    for file in cache_files
        filepath = joinpath(CACHE_DIR, file)
        size_bytes = stat(filepath).size
        size_gb = size_bytes / (1024^3)
        total_size += size_bytes
        
        push!(file_info, Dict(
            "file" => file,
            "size_gb" => size_gb,
            "modified" => stat(filepath).mtime
        ))
    end
    
    # Sort by modification time (newest first)
    sort!(file_info, by=f -> f["modified"], rev=true)
    
    return Dict(
        "exists" => true,
        "total_size_gb" => total_size / (1024^3),
        "num_files" => length(cache_files),
        "files" => file_info
    )
end

"""
Clean up old cache files if total size exceeds limit
"""
function cleanup_cache()
    info = get_cache_info()
    
    if !info["exists"] || info["total_size_gb"] <= MAX_CACHE_SIZE_GB
        return
    end
    
    files = info["files"]
    total_size = info["total_size_gb"]
    
    # Remove oldest files until we're under the limit
    files_to_remove = []
    for file_info in reverse(files)  # Start with oldest
        if total_size <= MAX_CACHE_SIZE_GB
            break
        end
        
        push!(files_to_remove, file_info["file"])
        total_size -= file_info["size_gb"]
    end
    
    # Remove the files
    for file in files_to_remove
        filepath = joinpath(CACHE_DIR, file)
        try
            rm(filepath)
            println("🗑️  Removed old cache file: $file")
        catch e
            @warn "Failed to remove cache file $file: $e"
        end
    end
end

"""
Clear all cached template grids
"""
function clear_cache()
    info = get_cache_info()
    
    if !info["exists"]
        return (0, 0.0)
    end
    
    removed_count = 0
    freed_gb = 0.0
    
    for file_info in info["files"]
        filepath = joinpath(CACHE_DIR, file_info["file"])
        try
            rm(filepath)
            removed_count += 1
            freed_gb += file_info["size_gb"]
        catch e
            @warn "Failed to remove cache file $(file_info["file"]): $e"
        end
    end
    
    return (removed_count, freed_gb)
end