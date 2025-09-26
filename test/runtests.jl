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
            @test (olig * olig2) == Olig("ACGTTGCA", "concat")
            @test description((olig * olig2)) == "concat"
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

            @testset "Display" begin
                olig = Olig("AGTC", "descr")
                olig_nodesc = Olig("AGTC", "")
                olig_long = Olig("A"^25, "")
                # Short display
                @test sprint(show, olig) == "Olig(\"AGTC\", len=4, desc=\"descr\")"
                @test sprint(show, olig_nodesc) == "Olig(\"AGTC\", len=4)"
                @test sprint(show, olig_long) == "Olig(\"AAAAAAAAAAAAAAAAA...\", len=25)"
                # Full display
                @test sprint(show, MIME"text/plain"(), olig) == "Olig\n  Sequence: AGTC\n  Length: 4\n  Description: \"descr\"\n"
                @test sprint(show, MIME"text/plain"(), olig_nodesc) == "Olig\n  Sequence: AGTC\n  Length: 4\n  Description: (none)\n"
                @test sprint(show, MIME"text/plain"(), olig_long) == "Olig\n  Sequence: AAAAAAAAAAAAAAAAAAAAAAAAA\n  Length: 25\n  Description: (none)\n"
            end
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

            @testset "Display" begin
                dolig = DegenerateOlig("AGNTC", "descr")
                dolig_nodesc = DegenerateOlig("AGNTC", "")
                dolig_long = DegenerateOlig("N"^25, "")
                # Short display
                @test sprint(show, dolig) == "DegenerateOlig(\"AGNTC\", len=5, n_deg=1, vars=4, desc=\"descr\")"
                @test sprint(show, dolig_nodesc) == "DegenerateOlig(\"AGNTC\", len=5, n_deg=1, vars=4)"
                @test sprint(show, dolig_long) == "DegenerateOlig(\"NNNNNNNNNNNNNNNNN...\", len=25, n_deg=25, vars=>10k)"
                # Full display
                @test sprint(show, MIME"text/plain"(), dolig) == "DegenerateOlig\n  Sequence: AGNTC\n  Length: 5\n  Degenerate positions: 1\n  Unique variants: 4\n  Description: \"descr\"\n"
                @test sprint(show, MIME"text/plain"(), dolig_nodesc) == "DegenerateOlig\n  Sequence: AGNTC\n  Length: 5\n  Degenerate positions: 1\n  Unique variants: 4\n  Description: (none)\n"
                @test sprint(show, MIME"text/plain"(), dolig_long) == "DegenerateOlig\n  Sequence: NNNNNNNNNNNNNNNNNNNNNNNNN\n  Length: 25\n  Degenerate positions: 25\n  Unique variants: 1125899906842624\n  Description: (none)\n"
            end
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
            @test subview isa AlignedOligs.OligView{Olig}
            @test n_deg_pos(subview) == 0
            @test n_unique_oligs(subview) == 1
            
            # View of degenerate view
            deg_subview = view2[2:3]
            @test String(deg_subview) == "GN"
            @test description(deg_subview) == "deg test"
            @test deg_subview isa AlignedOligs.OligView{DegenerateOlig}
            @test n_deg_pos(deg_subview) == 1
            @test n_unique_oligs(deg_subview) == 4

            @testset "Display" begin
                dolig_view = DegenerateOlig("AGNTC", "descr")[2:3]
                dolig_view_nodesc = DegenerateOlig("AGNTC", "")[2:4]
                dolig_view_long = DegenerateOlig("N"^25, "")[10:22]
                # Short display
                @test sprint(show, dolig_view) == "OligView(\"GN\", len=2, range=2:3, desc=\"descr\")"
                @test sprint(show, dolig_view_nodesc) == "OligView(\"GNT\", len=3, range=2:4)"
                @test sprint(show, dolig_view_long) == "OligView(\"NNNNNNNNNNNNN\", len=13, range=10:22)"
                # Full display
                @test sprint(show, MIME"text/plain"(), dolig_view) == "OligView{DegenerateOlig}\n  Viewed sequence: GN\n  Length: 2\n  Range: 2:3\n  Parent description: descr\n"
                @test sprint(show, MIME"text/plain"(), dolig_view_nodesc) == "OligView{DegenerateOlig}\n  Viewed sequence: GNT\n  Length: 3\n  Range: 2:4\n  Parent description: \n"
                @test sprint(show, MIME"text/plain"(), dolig_view_long) == "OligView{DegenerateOlig}\n  Viewed sequence: NNNNNNNNNNNNN\n  Length: 13\n  Range: 10:22\n  Parent description: \n"
            end
        end
        
        @testset "GappedOlig" begin
            # Construction tests
            olig = Olig("ACGTACGT", "gapped test")
            gaps = [3 => 2, 7 => 1]  # positions 3-4 gap of 2, position 7 gap of 1
            go = GappedOlig(olig, gaps, "gapped test")
            @test go == GappedOlig("AC--GT-ACGT")
            @test go isa GappedOlig
            @test go.olig == olig
            @test go.gaps == gaps
            @test go.total_length == 8 + 2 + 1 == 11

            # Invalid constructions
            @test_throws ArgumentError GappedOlig(olig, [3 => 0])  # zero length gap
            @test_throws ArgumentError GappedOlig(olig, [3 => 2, 4 => 1])  # overlapping
            @test_throws ArgumentError GappedOlig(olig, [10 => 1])  # exceeds length
            @test_throws ArgumentError GappedOlig(olig, [0 => 1])  # invalid start (assuming added check)
            @test_throws ArgumentError GappedOlig(olig, [-1 => 1])  # negative start
            @test GappedOlig(olig, Pair{Int,Int}[], "no gaps") isa GappedOlig  # empty gaps allowed

            # String representation
            @test String(go) == "AC--GT-ACGT"
            @test String(GappedOlig(olig, Pair{Int,Int}[])) == "ACGTACGT"

            # Length and isempty
            @test length(go) == 11
            @test isempty(go) == false
            empty_go = GappedOlig(Olig(""), Pair{Int,Int}[], "empty")
            @test isempty(empty_go) == true
            @test length(empty_go) == 0
            @test String(empty_go) == ""
            gap_only = GappedOlig(Olig(""), [1 => 3], "gap only")
            @test length(gap_only) == 3
            @test String(gap_only) == "---"

            # Indexing and char_at
            @test go[1] == 'A'
            @test go[2] == 'C'
            @test go[3] == '-'
            @test go[4] == '-'
            @test go[5] == 'G'
            @test go[6] == 'T'
            @test go[7] == '-'
            @test go[8] == 'A'
            @test go[9] == 'C'
            @test go[10] == 'G'
            @test go[11] == 'T'
            @test_throws BoundsError go[0]
            @test_throws BoundsError go[12]

            # Iteration
            @test collect(go) == ['A','C','-','-','G','T','-','A','C','G','T']
            @test first(go) == 'A'
            @test last(go) == 'T'

            # Slicing (as OligView)
            @test String(go[3:7]) == "--GT-"
            @test String(go[1:5]) == "AC--G"
            @test String(go[8:end]) == "ACGT"

            # Equality
            go2 = GappedOlig(Olig("ACGTACGT"), [3=>2,7=>1])
            @test go == go2
            @test go == "AC--GT-ACGT"
            @test "AC--GT-ACGT" == go
            @test go != GappedOlig(Olig("ACGTACGT"), [3=>1])

            # Description
            @test description(go) == "gapped test"
            @test description(GappedOlig(olig, Pair{Int,Int}[])) == "gapped test"
            @test description(GappedOlig(olig, Pair{Int,Int}[], "new descr")) == "new descr"

            # hasgaps
            @test hasgaps(go) == true
            @test hasgaps(GappedOlig(olig, Pair{Int,Int}[])) == false

            # Properties
            @test n_deg_pos(go) == 0
            @test n_unique_oligs(go) == 1

            # revcomp
            rev_go = SeqFold.revcomp(go)
            @test rev_go isa GappedOlig
            @test rev_go.gaps == [5=>1, 8=>2]  # reversed positions
            @test String(rev_go) == "ACGT-AC--GT"  # verified revcomp with gaps adjusted

            # complement
            comp_go = SeqFold.complement(go)
            @test comp_go isa GappedOlig
            @test comp_go.gaps == go.gaps  # gaps unchanged
            @test String(comp_go) == "TG--CA-TGCA"  # complement of underlying, gaps preserved

            # gc_content
            @test SeqFold.gc_content(go) ≈ SeqFold.gc_content(olig)  # ignores gaps
            @test isnan(SeqFold.gc_content(GappedOlig(Olig(""), [1=>1])))  # empty underlying

            # Folding and energy errors if gaps present
            @test_throws ErrorException SeqFold.fold(go)
            @test_throws ErrorException SeqFold.dg(go)
            @test_throws ErrorException SeqFold.dg_cache(go)
            @test_throws ErrorException SeqFold.tm(go, go)
            @test_throws ErrorException SeqFold.tm_cache(go, go)

            # Non-gapped cases don't throw
            no_gap_go = GappedOlig(olig, Pair{Int,Int}[])
            @test SeqFold.fold(no_gap_go) == SeqFold.fold(olig.seq)  # assuming fold takes String
            @test SeqFold.dg(no_gap_go) == SeqFold.dg(olig.seq)

            # OligView on GappedOlig
            go_view = go[2:6]
            @test go_view isa AlignedOligs.OligView{GappedOlig}
            @test String(go_view) == "C--GT"
            @test length(go_view) == 5
            @test collect(go_view) == ['C','-','-','G','T']
            @test description(go_view) == "gapped test"
            @test n_deg_pos(go_view) == 0
            @test n_unique_oligs(go_view) == 1

            # Subview
            sub_go_view = go_view[2:4]
            @test String(sub_go_view) == "--G"
            @test length(sub_go_view) == 3

            # NonDegenIterator
            @test length(nondegens(go)) == 1
            @test collect(nondegens(go)) == [go]

            # Empty GappedOlig
            @test collect(empty_go) == Char[]
            @test n_deg_pos(empty_go) == 0
            @test n_unique_oligs(empty_go) == 1


            @testset "Display" begin
                golig = GappedOlig("AG-T-C", "descr")
                golig_nodesc = GappedOlig("AG-T-C", "")
                golig_long = GappedOlig("A-"^10, "")
                # Short display
                @test sprint(show, golig) == "GappedOlig(\"AG-T-C\", len=6, gaps=2, desc=\"descr\")"
                @test sprint(show, golig_nodesc) == "GappedOlig(\"AG-T-C\", len=6, gaps=2)"
                @test sprint(show, golig_long) == "GappedOlig(\"A-A-A-A-A-A-A-A-A-A-\", len=20, gaps=10)"
                # Full display
                @test sprint(show, MIME"text/plain"(), golig) == "GappedOlig\n  Gapped sequence: AG-T-C\n  Length (with gaps): 6\n  Underlying Olig: AGTC\n  Gaps: [3 => 1, 5 => 1]\n  Description: \"descr\"\n"
                @test sprint(show, MIME"text/plain"(), golig_nodesc) == "GappedOlig\n  Gapped sequence: AG-T-C\n  Length (with gaps): 6\n  Underlying Olig: AGTC\n  Gaps: [3 => 1, 5 => 1]\n  Description: (none)\n"
                @test sprint(show, MIME"text/plain"(), golig_long) == "GappedOlig\n  Gapped sequence: A-A-A-A-A-A-A-A-A-A-\n  Length (with gaps): 20\n  Underlying Olig: AAAAAAAAAA\n  Gaps: [2 => 1, 4 => 1, 6 => 1, 8 => 1, 10 => 1, 12 => 1, 14 => 1, 16 => 1, 18 => 1, 20 => 1]\n  Description: (none)\n"
            end
        end
    end
end
