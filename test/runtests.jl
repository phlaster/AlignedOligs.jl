using Test
using Aqua
using JET
using Random

using AlignedOligs
using AlignedOligs.Oligs
using AlignedOligs.Alignments
using AlignedOligs.Primers

Random.seed!(42)

@testset verbose=true failfast=true "AlignedOligs.jl"  begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(AlignedOligs)
    end
    @testset "Code linting (JET.jl)" begin
        JET.test_package(AlignedOligs; target_defined_modules = true)
    end
    
    @testset "Oligs" include("test_oligs.jl")

    # @testset "Alignments" include("test_alignments.jl")
    # @testset "Primers" include("test_primers.jl")
    # @testset "SeqFold methods" include("test_seqfold.jl")


end