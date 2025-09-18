using AlignedOligs
using Test
using Aqua
using JET

@testset "AlignedOligs.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(AlignedOligs)
    end
    @testset "Code linting (JET.jl)" begin
        JET.test_package(AlignedOligs; target_defined_modules = true)
    end
    # Write your tests here.
end
