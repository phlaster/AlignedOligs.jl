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
    
    @testset "Oligs" begin
    
        @testset "Olig" begin
            # Construction tests
            @test Olig("ACGT") isa Olig
            @test Olig("ACGT", "test description") isa Olig
            @test Olig("ACGT", "") isa Olig
            @test Olig("ACGT", nothing) isa Olig
            @test Olig(['A', 'C', 'G', 'T']) isa Olig
            @test Olig(['A', 'C', 'G', 'T'], "test description") isa Olig
            @test Olig("ACGT", 123) isa Olig
            @test_throws ErrorException Olig("ACGTX")
            @test Olig() == Olig("")
            @test Olig("", "empty") == Olig("", "empty")
            
            # String interface tests
            olig = Olig("ACGT", "test")
            @test String(olig) == "ACGT"
            @test length(olig) == 4
            @test isempty(Olig("")) == true
            @test isempty(olig) == false
            @test collect(olig) == ['A', 'C', 'G', 'T']
            @test olig[1] == 'A'
            @test olig[2:3] == "CG"
            @test occursin("CG", olig)
            @test first(olig) == 'A'
            @test last(olig) == 'T'
            @test olig[4:-1:1] == "TGCA"
            
            # Description tests
            @test description(olig) == "test"
            @test description(Olig("ACGT", "")) == ""
            
            # Conversion tests
            deg_olig = DegenerateOlig("ACGT")
            @test Olig(deg_olig) isa Olig
            @test_throws InexactError Olig(DegenerateOlig("ACGN"))
            @test convert(Olig, deg_olig) isa Olig
            
            # Concatenation tests
            olig2 = Olig("TGCA")
            @test (olig * olig2) == Olig("ACGTTGCA", "concatenated")
            @test description((olig * olig2)) == "concatenated"
            @test (Olig("A") * Olig("C") * Olig("G") * Olig("T")) == olig
            
            # Equality tests
            @test Olig("ACGT") == Olig("ACGT")
            @test Olig("ACGT") != Olig("TGCA")
            @test Olig("ACGT") == "ACGT"
            @test "ACGT" == Olig("ACGT")
            @test Olig("ACGT") != "TGCA"
            
            # Case handling
            @test Olig("acgt") == Olig("ACGT")
            
            # Empty sequence tests
            empty_olig = Olig("", "empty")
            @test isempty(empty_olig) == true
            @test length(empty_olig) == 0
            @test String(empty_olig) == ""
            @test_throws BoundsError empty_olig[1]
        end
        
        @testset "DegenerateOlig" begin
            # Construction tests
            @test DegenerateOlig("ACGT") isa DegenerateOlig
            @test DegenerateOlig("ACGN") isa DegenerateOlig
            @test DegenerateOlig("ACGN", "test description") isa DegenerateOlig
            @test DegenerateOlig("ACGN", "") isa DegenerateOlig
            @test DegenerateOlig("ACGN", nothing) isa DegenerateOlig
            @test DegenerateOlig(['A', 'C', 'G', 'N']) isa DegenerateOlig
            @test DegenerateOlig(['A', 'C', 'G', 'N'], "test description") isa DegenerateOlig
            @test DegenerateOlig("ACGT", 123) isa DegenerateOlig
            @test_throws ErrorException DegenerateOlig("ACGTX")
            @test DegenerateOlig() == DegenerateOlig("")
            @test DegenerateOlig("", "empty") == DegenerateOlig("", 0, 1, "empty")
            
            # String interface tests
            deg_olig = DegenerateOlig("ACGN", "test")
            @test String(deg_olig) == "ACGN"
            @test length(deg_olig) == 4
            @test isempty(DegenerateOlig(""))
            @test !isempty(deg_olig)
            @test collect(deg_olig) == ['A', 'C', 'G', 'N']
            @test deg_olig[1] == 'A'
            @test deg_olig[2:3] == "CG"
            @test occursin("CG", deg_olig)
            @test first(deg_olig) == 'A'
            @test last(deg_olig) == 'N'
            @test deg_olig[4:-1:1] == "NGCA"
            
            # Properties tests
            @test n_deg_pos(deg_olig) == 1
            @test n_unique_oligs(deg_olig) == 4
            @test n_deg_pos(DegenerateOlig("ACGT")) == 0
            @test n_unique_oligs(DegenerateOlig("ACGT")) == 1
            @test n_deg_pos(DegenerateOlig("NNNN")) == 4
            @test n_unique_oligs(DegenerateOlig("NNNN")) == 256
            @test n_deg_pos(DegenerateOlig("RY")) == 2
            @test n_unique_oligs(DegenerateOlig("RY")) == 4
            @test n_deg_pos(DegenerateOlig("BVDH")) == 4
            @test n_unique_oligs(DegenerateOlig("BVDH")) == 81
            
            # Description tests
            @test description(deg_olig) == "test"
            @test description(DegenerateOlig("ACGN", "")) == ""
            
            # Conversion tests
            olig = Olig("ACGT")
            @test DegenerateOlig(olig) isa DegenerateOlig
            @test convert(DegenerateOlig, olig) isa DegenerateOlig
            @test String(DegenerateOlig(olig)) == "ACGT"
            @test n_deg_pos(DegenerateOlig(olig)) == 0
            @test n_unique_oligs(DegenerateOlig(olig)) == 1
            
            # Concatenation tests
            deg_olig2 = DegenerateOlig("NNGT")
            @test (deg_olig * deg_olig2) isa DegenerateOlig
            @test String((deg_olig * deg_olig2)) == "ACGNNNGT"
            @test n_deg_pos((deg_olig * deg_olig2)) == 3
            @test n_unique_oligs((deg_olig * deg_olig2)) == 4^3
            
            @test (olig * deg_olig) isa DegenerateOlig
            @test String((olig * deg_olig)) == "ACGTACGN"
            @test n_deg_pos((olig * deg_olig)) == 1
            @test n_unique_oligs((olig * deg_olig)) == 4
            
            @test (deg_olig * olig) isa DegenerateOlig
            @test String((deg_olig * olig)) == "ACGNACGT"
            @test n_deg_pos((deg_olig * olig)) == 1
            @test n_unique_oligs((deg_olig * olig)) == 4
            
            # NonDegenIterator tests
            @test length(nondegens(deg_olig)) == 4
            variants = collect(nondegens(deg_olig))
            @test length(variants) == 4
            @test Set(String.(variants)) == Set(["ACGA", "ACGC", "ACGG", "ACGT"])
            
            # Complex degenerate sequence
            complex_deg = DegenerateOlig("RYSWKMBDHVN")
            @test n_deg_pos(complex_deg) == 11
            @test n_unique_oligs(complex_deg) == 20736
            
            # NonDegenIterator for complex sequence
            @test length(nondegens(complex_deg)) == 20736
            complex_variants = collect(Iterators.take(nondegens(complex_deg), 5))
            @test all(x -> length(x) == 11, complex_variants)
            
            # Empty sequence
            empty_deg = DegenerateOlig("", "empty")
            @test isempty(empty_deg) == true
            @test length(empty_deg) == 0
            @test String(empty_deg) == ""
            @test n_deg_pos(empty_deg) == 0
            @test n_unique_oligs(empty_deg) == 1
            @test length(nondegens(empty_deg)) == 1
        end
        
        @testset "OligView" begin
            # Basic view creation
            olig = Olig("ACGTACGT", "test")
            deg_olig = DegenerateOlig("ACGNACGN", "deg test")
            
            view1 = olig[2:5]
            @test view1 isa AlignedOligs.OligView{Olig}
            @test String(view1) == "CGTA"
            @test length(view1) == 4
            @test collect(view1) == ['C', 'G', 'T', 'A']
            @test description(view1) == "test"
            
            view2 = deg_olig[2:5]
            @test view2 isa AlignedOligs.OligView{DegenerateOlig}
            @test String(view2) == "CGNA"
            @test length(view2) == 4
            @test collect(view2) == ['C', 'G', 'N', 'A']
            @test description(view2) == "deg test"
            
            # Edge cases
            @test isempty(olig[1:0])
            @test_throws BoundsError olig[10:11]
            
            # String interface
            @test view1[1] == 'C'
            @test view1[2:3] == "GT"
            @test view1[end] == 'A'
            @test view1[3:-1:1] == "TGC"
            @test occursin("GT", view1)
            
            # Properties
            @test n_deg_pos(view1) == 0
            @test n_unique_oligs(view1) == 1
            @test n_deg_pos(view2) == 1
            @test n_unique_oligs(view2) == 4
            
            # Concatenation
            @test (view1 * Olig("TG")) isa Olig
            @test String(view1 * Olig("TG")) == "CGTATG"
            
            @test (view2 * Olig("TG")) isa DegenerateOlig
            @test String(view2 * Olig("TG")) == "CGNATG"
            @test n_deg_pos(view2 * Olig("TG")) == 1
            @test n_unique_oligs(view2 * Olig("TG")) == 4
            
            @test (Olig("AT") * view1) isa Olig
            @test String(Olig("AT") * view1) == "ATCGTA"
            
            @test (Olig("AT") * view2) isa DegenerateOlig
            @test String(Olig("AT") * view2) == "ATCGNA"
            @test n_deg_pos(Olig("AT") * view2) == 1
            @test n_unique_oligs(Olig("AT") * view2) == 4
            
            # Conversion
            @test convert(Olig, view1) isa Olig
            @test String(convert(Olig, view1)) == "CGTA"
            @test_throws InexactError convert(Olig, view2)
            
            # NonDegenIterator for views
            @test length(nondegens(view2)) == 4
            view_variants = collect(nondegens(view2))
            @test length(view_variants) == 4
            @test Set(String.(view_variants)) == Set(["CGAA", "CGCA", "CGGA", "CGTA"])
            
            # Empty view
            empty_view = olig[1:0]
            @test isempty(empty_view)
            @test length(empty_view) == 0
            @test collect(empty_view) == Char[]
            @test n_deg_pos(empty_view) == 0
            @test n_unique_oligs(empty_view) == 1
            
            # View of a view
            subview = view1[2:3]
            @test String(subview) == "GT"
            @test description(subview) == "test"
            @test subview isa AlignedOligs.OligView{AlignedOligs.OligView{Olig}}
            @test n_deg_pos(subview) == 0
            @test n_unique_oligs(subview) == 1
            
            # View of degenerate view
            deg_subview = view2[2:3]
            @test String(deg_subview) == "GN"
            @test description(deg_subview) == "deg test"
            @test deg_subview isa AlignedOligs.OligView{AlignedOligs.OligView{DegenerateOlig}}
            @test n_deg_pos(deg_subview) == 1
            @test n_unique_oligs(deg_subview) == 4
        end
    end
end
