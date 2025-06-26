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
end
