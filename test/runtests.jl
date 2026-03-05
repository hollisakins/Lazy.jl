using Lazy
using Test

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
end
