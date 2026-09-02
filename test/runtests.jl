using Lazy
using Test
using Trapz
using HDF5
using FITSIO

@testset "Lazy.jl" begin
    @testset "Core functionality" begin
        # Test LazyError type
        err = LazyError("test message")
        @test err.msg == "test message"
        @test isa(err, Exception)
        
        # Test that template/filter directories exist (should pass normally)
        @test isdir(Lazy.templatepath)
        @test isdir(Lazy.filterpath)
    end
    
    @testset "CLI argument parsing" begin
        # Test help output doesn't crash
        @test main(String[]) == 0
        
        # Test invalid arguments throw LazyError
        @test_throws LazyError main(String["invalid"])
        @test_throws LazyError main(String["fit", "invalid"])
        @test_throws LazyError main(String["fit", "-p"])  # Missing param file
    end

    @testset "CGM damping wing (Asada+24)" begin
        A, a, c = 3.5918, 1.8414, 18.001
        @test Lazy.cgm_sigmoid(6.0, A, a, c) ≈ A / 2.0 + c
        @test Lazy.cgm_sigmoid(20.0, A, a, c) ≈ A + c atol=1e-4
        @test Lazy.cgm_sigmoid(0.0, A, a, c) ≈ c atol=1e-4

        nu_lya = 2.46607e15
        @test Lazy.sigma_a(nu_lya * 0.99) > 0.0
        @test isfinite(Lazy.sigma_a(nu_lya * 0.99))
        @test isfinite(Lazy.sigma_a(nu_lya * 1.01))
        @test Lazy.sigma_a(nu_lya * 0.99) > Lazy.sigma_a(nu_lya * 0.5)

        templwav = collect(range(900.0, 2000.0, length=500))
        t_low = ones(length(templwav))
        Lazy.apply_cgm_damping_wing!(t_low, templwav, 5.0, A, a, c)
        @test all(t_low .== 1.0)

        t_high = ones(length(templwav))
        Lazy.apply_cgm_damping_wing!(t_high, templwav, 8.0, A, a, c)
        @test all(t_high .<= 1.0)
        @test any(t_high .< 1.0)

        idx_lya = argmin(abs.(templwav .- 1216.0))
        idx_far = argmin(abs.(templwav .- 1800.0))
        @test t_high[idx_lya] < t_high[idx_far]
    end

    @testset "Spectroscopic redshift fixing (use_zspec)" begin
        nz_test, nband_test, ntempl_test = 5, 2, 1
        zgrid_test = collect(range(0.0, 2.0, length=nz_test))
        templgrid_test = ones(Float32, nband_test, ntempl_test, nz_test)
        template_error_test = zeros(nz_test, nband_test)
        fnu_test = [1.0, 1.0]
        efnu_test = [0.1, 0.1]

        zbest_free, _, _, chi2_row_free = Lazy.fit_single_object(
            1, fnu_test, efnu_test, templgrid_test, template_error_test,
            zgrid_test, 2, nband_test, ntempl_test, nz_test)
        @test zbest_free >= 0
        @test count(chi2_row_free .< 1e18) == nz_test

        zbest_fix, _, _, chi2_row_fix = Lazy.fit_single_object(
            1, fnu_test, efnu_test, templgrid_test, template_error_test,
            zgrid_test, 2, nband_test, ntempl_test, nz_test;
            z_fix_idx=3)
        @test zbest_fix == zgrid_test[3]
        @test chi2_row_fix[3] < 1e18
        @test all(chi2_row_fix[i] == Float32(1e18) for i in [1,2,4,5])
    end

    @testset "P(z) from high-chi2 grids (underflow regression)" begin
        zgrid = collect(0.01:0.01:20.0)
        nz = length(zgrid)

        # Poor fit: min chi2 = 178 used to underflow exp(-chi2/2) and
        # overflow Float32(1/total), producing all-NaN/Inf P(z)
        chi2_poor = Float32.(178.0 .+ 2000.0 .* (zgrid .- 7.0) .^ 2 ./ (1 .+ zgrid))
        pz = Lazy.chi2grid_to_pz(reshape(chi2_poor, nz, 1), zgrid)[:, 1]
        @test all(isfinite, pz)
        @test trapz(zgrid, pz) ≈ 1.0 atol=1e-6
        @test argmax(pz) == argmin(chi2_poor)
        @test all(isfinite, Float32.(pz))  # survives the Float32 output conversion

        # Even chi2 in the thousands (underflows Float64 without the shift)
        chi2_awful = Float32.(5000.0 .+ 100.0 .* (zgrid .- 10.0) .^ 2)
        pz_awful = Lazy.chi2grid_to_pz(reshape(chi2_awful, nz, 1), zgrid)[:, 1]
        @test all(isfinite, pz_awful)
        @test trapz(zgrid, pz_awful) ≈ 1.0 atol=1e-6

        # Shift-invariance: a good fit gives the same P(z) as the unshifted formula
        chi2_good = Float32.(3.0 .+ 50.0 .* (zgrid .- 2.0) .^ 2)
        pz_ref = exp.(-0.5 .* Float64.(chi2_good))
        pz_ref ./= trapz(zgrid, pz_ref)
        pz_good = Lazy.chi2grid_to_pz(reshape(chi2_good, nz, 1), zgrid)[:, 1]
        @test pz_good ≈ pz_ref rtol=1e-6

        # Negative chi2 is an invalid-fit sentinel: P(z) = 0 there
        chi2_sent = copy(chi2_good)
        chi2_sent[1:100] .= -1.0f0
        pz_sent = Lazy.chi2grid_to_pz(reshape(chi2_sent, nz, 1), zgrid)[:, 1]
        @test all(pz_sent[1:100] .== 0.0)
        @test trapz(zgrid, pz_sent) ≈ 1.0 atol=1e-6

        # A column with no valid chi2 stays all-zero
        pz_none = Lazy.chi2grid_to_pz(fill(-1.0f0, nz, 2), zgrid)
        @test all(pz_none .== 0.0)

        # Inf/NaN chi2 (e.g. from old work files) get P(z) = 0, not NaN
        chi2_dirty = copy(chi2_good)
        chi2_dirty[1] = Inf32
        chi2_dirty[2] = NaN32
        pz_dirty = Lazy.chi2grid_to_pz(reshape(chi2_dirty, nz, 1), zgrid)[:, 1]
        @test pz_dirty[1] == 0.0 && pz_dirty[2] == 0.0
        @test all(isfinite, pz_dirty)
        @test trapz(zgrid, pz_dirty) ≈ 1.0 atol=1e-6
        pz_allinf = Lazy.chi2grid_to_pz(fill(Inf32, nz, 1), zgrid)
        @test all(pz_allinf .== 0.0)

        # z_spec-mode filler (1e18) and poor-fit marker (1e10) are sentinels,
        # not real chi2: they must not become the shift minimum. A rejected
        # fixed-z fit (filler everywhere, -1 at the fixed bin) has no P(z).
        chi2_zfix_rejected = fill(Lazy.CHI2_ZFIX_FILLER, nz)
        chi2_zfix_rejected[500] = -1.0f0
        pz_zfix = Lazy.chi2grid_to_pz(reshape(chi2_zfix_rejected, nz, 1), zgrid)[:, 1]
        @test all(pz_zfix .== 0.0)

        # A successful fixed-z fit yields a delta function at the fitted bin
        chi2_zfix_ok = fill(Lazy.CHI2_ZFIX_FILLER, nz)
        chi2_zfix_ok[500] = 250.0f0
        pz_zfix_ok = Lazy.chi2grid_to_pz(reshape(chi2_zfix_ok, nz, 1), zgrid)[:, 1]
        @test all(isfinite, pz_zfix_ok)
        @test pz_zfix_ok[500] > 0
        @test all(pz_zfix_ok[[1:499; 501:nz]] .== 0.0)

        # All-poor-fit column (marker at every z) has no P(z)
        pz_poor = Lazy.chi2grid_to_pz(fill(Lazy.CHI2_POOR_FIT, nz, 1), zgrid)
        @test all(pz_poor .== 0.0)
    end

    @testset "HDF5 -> FITS export (pure-Julia CFITSIO writer)" begin
        # Build a tiny, fully-populated work file and round-trip it through
        # convert_hdf5_to_fits, checking the SUMMARY / PZ / TEMPL layout that
        # downstream readers (astropy etc.) rely on.
        mktempdir() do dir
            work_file = joinpath(dir, "work.h5")
            fits_file = joinpath(dir, "out.fits")

            nobj, ntempl = 3, 2
            zgrid = collect(0.0:0.5:2.0)          # nz = 5, max z = 2
            nz = length(zgrid)
            bands = ["f150w", "f277w"]
            nband = length(bands)
            templates = [joinpath(dir, "templ_a.dat"), joinpath(dir, "templ_b.fits")]
            param = Dict{String,Any}(
                "io" => Dict{String,Any}("output_pz" => true, "output_templates" => true),
                "fitting" => Dict{String,Any}("use_zspec" => true),
            )

            Lazy.create_hdf5_work_file(work_file, param, nobj, nz, nband, ntempl,
                                       bands, zgrid, templates)

            ids = Int32[11, 22, 33]
            zbest = [0.5, 1.0, 1.5]
            coeffs = [1.0 2.0; 3.0 4.0; 5.0 6.0]        # (nobj, ntempl)
            chi2grid = Float32[i + 10j for i in 1:nz, j in 1:nobj]  # (nz, nobj)
            pz = Float32.(Lazy.chi2grid_to_pz(chi2grid, zgrid))
            h5open(work_file, "r+") do f
                f["results/ID"][:] = ids
                f["results/zbest"][:] = zbest
                f["results/chi2best"][:] = [1.0, 2.0, 3.0]
                f["results/z_spec"][:] = [0.4, 1.1, -1.0]
                for k in ["z_l95", "z_l68", "z_med", "z_u68", "z_u95"]
                    f["results/$k"][:] = zbest
                end
                f["results/coeffs"][:, :] = coeffs
                f["results/pz_gt1"][:] = Float32[0.1, 0.2, 0.3]
                f["results/pz_gt2"][:] = Float32[0.0, 0.0, 0.0]
                f["results/Sz"][:] = [0.9, 0.8, 0.7]
                for bc in 0:2
                    f["results/Pz$bc"][:] = Float32[bc, bc, bc]
                end
                f["results/Pz_cen"][:] = Float32[0.5, 0.6, 0.7]
                f["results/Pz_zgtrzb2"][:] = Float32[0.05, 0.06, 0.07]
                f["photometry/f150w"][:] = [10.0, 20.0, 30.0]
                f["photometry/f277w"][:] = [11.0, 21.0, 31.0]
                f["pz/chi2grid"][:, :] = chi2grid
                f["pz/pz"][:, :] = pz
            end

            wav = collect(1000.0:500.0:3000.0)        # nwav = 5
            nwav = length(wav)
            templgrid = zeros(Float32, nwav, ntempl, nz)
            for t in 1:ntempl, iz in 1:nz
                templgrid[:, t, iz] .= Float32.(100t .+ iz .+ (1:nwav) ./ 10)
            end
            Lazy.write_template_data(work_file, templgrid, wav, zgrid, templates)
            Lazy.finalize_hdf5_work_file(work_file)

            Lazy.convert_hdf5_to_fits(work_file, fits_file; chunk_size=2)
            @test isfile(fits_file)

            FITS(fits_file, "r") do ff
                @test length(ff) == 4  # primary + SUMMARY + PZ + TEMPL
                @test read_key(ff[2], "EXTNAME")[1] == "SUMMARY"
                @test read_key(ff[3], "EXTNAME")[1] == "PZ"
                @test read_key(ff[4], "EXTNAME")[1] == "TEMPL"

                # ── SUMMARY ──
                summ = ff[2]
                cols = FITSIO.colnames(summ)
                @test cols[1:4] == ["ID", "z_best", "chi2", "z_spec"]
                for c in ["z_l95", "z_l68", "z_med", "z_u68", "z_u95", "Pz_gt1", "Pz_gt2",
                          "Sz", "Pz0", "Pz1", "Pz2", "Pz_cen", "Pz_zgtrzb2",
                          "f150w", "f277w", "coeffs"]
                    @test c in cols
                end
                @test read(summ, "ID") == Int64.(ids)
                @test read(summ, "z_best") == Float32.(zbest)
                @test read(summ, "z_spec") == Float32[0.4, 1.1, -1.0]
                @test read(summ, "Pz_gt1") == Float32[0.1, 0.2, 0.3]
                @test read(summ, "Pz2") == Float32[2, 2, 2]
                @test read(summ, "f277w") == Float32[11, 21, 31]
                # coeffs is an ntempl-wide vector column: FITSIO returns (ntempl, nobj)
                @test read(summ, "coeffs") == Float32.(permutedims(coeffs))
                iband = findfirst(==("f150w"), cols)
                @test read_key(summ, "TUNIT$iband")[1] == "fnu"
                icoef = findfirst(==("coeffs"), cols)
                @test read_key(summ, "TFORM$icoef")[1] == "$(ntempl)E"

                # ── PZ: sentinel row (ID=-1, zgrid) followed by one row per object ──
                pzh = ff[3]
                @test FITSIO.colnames(pzh) == ["ID", "Pz"]
                @test read(pzh, "ID") == vcat(Int64[-1], Int64.(ids))
                pz_out = read(pzh, "Pz")                # (nz, nobj+1)
                @test size(pz_out) == (nz, nobj + 1)
                @test pz_out[:, 1] == Float32.(zgrid)
                @test pz_out[:, 2:end] == pz

                # ── TEMPL: sentinel row (z=-1, wavelength) followed by one row per z ──
                th = ff[4]
                @test FITSIO.colnames(th) == ["z", "templ_a", "templ_b"]
                @test read(th, "z") == vcat(Float32[-1], Float32.(zgrid))
                ta = read(th, "templ_a")               # (nwav, nz+1)
                @test size(ta) == (nwav, nz + 1)
                @test ta[:, 1] == Float32.(wav)
                @test ta[:, 2:end] == templgrid[:, 1, :]
                tb = read(th, "templ_b")
                @test tb[:, 2:end] == templgrid[:, 2, :]
            end
        end
    end
end
