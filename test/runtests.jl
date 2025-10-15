using AlignedOligs
using Test
using Aqua
using JET

@testset "AlignedOligs.jl" begin
    # @testset "Code quality (Aqua.jl)" begin
    #     Aqua.test_all(AlignedOligs)
    # end
    # @testset "Code linting (JET.jl)" begin
    #     JET.test_package(AlignedOligs; target_defined_modules = true)
    # end
    
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
            @test view1 isa OligView{Olig}
            @test String(view1) == "CGTA"
            @test length(view1) == 4
            @test collect(view1) == ['C', 'G', 'T', 'A']
            @test description(view1) == "test"
            
            view2 = deg_olig[2:5]
            @test view2 isa OligView{DegenerateOlig}
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
            @test subview isa OligView{Olig}
            @test n_deg_pos(subview) == 0
            @test n_unique_oligs(subview) == 1
            
            # View of degenerate view
            deg_subview = view2[2:3]
            @test String(deg_subview) == "GN"
            @test description(deg_subview) == "deg test"
            @test deg_subview isa OligView{DegenerateOlig}
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
            # Construction tests with Olig parent
            seq = "ACGTACGT"
            gapped_seq = "AC--GT-ACGT"
            olig = Olig(seq, "test")
            go = GappedOlig(olig, [3=>2, 5=>1])
            @test go isa GappedOlig{Olig}
            @test String(go) == "AC--GT-ACGT"
            @test length(go) == 11
            @test parent(go) === olig
            @test go.gaps == [3=>2, 5=>1]
            @test description(go) == description(olig)
            @test hasgaps(go)
            @test n_deg_pos(go) == 0
            @test n_unique_oligs(go) == 1
            @test collect(go) == collect(gapped_seq)
            
            # Construction from string
            go_str = GappedOlig(gapped_seq, "test")
            @test go_str isa GappedOlig{Olig}
            @test String(go_str) == "AC--GT-ACGT"
            @test go_str.parent == Olig(seq, "test")
            @test go_str.gaps == [3=>2, 5=>1]
            @test length(go_str) == 11
            
            # Construction with DegenerateOlig parent
            deg_seq = "ACGNRT"
            gapped_deg_seq = "AC--GN-RT"
            deg_olig = DegenerateOlig(deg_seq, "deg test")
            go_deg = GappedOlig(deg_olig, [3=>2, 5=>1])
            @test go_deg isa GappedOlig{DegenerateOlig}
            @test String(go_deg) == gapped_deg_seq
            @test length(go_deg) == 9
            @test parent(go_deg) === deg_olig
            @test go_deg.gaps == [3=>2, 5=>1]
            @test n_deg_pos(go_deg) == 2
            @test n_unique_oligs(go_deg) == 8
            @test collect(go_deg) == collect(gapped_deg_seq)
            
            # Construction from degenerate string
            go_deg_str = GappedOlig(gapped_deg_seq, "deg test")
            @test go_deg_str isa GappedOlig{DegenerateOlig}
            @test String(go_deg_str) == gapped_deg_seq
            @test go_deg_str.parent == DegenerateOlig(deg_seq, "deg test")
            @test go_deg_str.gaps == [3=>2, 5=>1]
            
            # Invalid gap positions
            @test_throws ArgumentError GappedOlig(olig, [0=>1]) # Start < 1
            @test_throws ArgumentError GappedOlig(olig, [10=>1]) # Start > parent len + 1
            @test_throws ArgumentError GappedOlig(olig, [3=>0]) # Non-positive gap length
            @test_throws ArgumentError GappedOlig(olig, [3=>2, 3=>1]) # Overlapping gaps
            
            # Edge cases: empty and no gaps
            empty_go = GappedOlig(Olig(""), Pair{Int64, Int64}[])
            @test String(empty_go) == ""
            @test length(empty_go) == 0
            @test isempty(empty_go) == true
            @test hasgaps(empty_go) == false
            no_gap_go = GappedOlig(olig, Pair{Int64, Int64}[])
            @test String(no_gap_go) == "ACGTACGT"
            @test length(no_gap_go) == 8
            @test hasgaps(no_gap_go) == false
            
            # Indexing tests
            @test go[1] == 'A'
            @test go[3] == '-'
            @test go[5] == 'G'
            @test go[11] == 'T'
            @test_throws BoundsError go[12]
            @test go[1:3] isa OligView{GappedOlig{Olig}}
            @test String(go[1:3]) == "AC-"
            @test String(go[3:4]) == "--"
            @test String(go[5:7]) == "GT-"
            @test String(go[8:11]) == "ACGT"
            @test String(go[11:11]) == "T"
            
            # Slicing with DegenerateOlig
            go_view = go_deg[2:7]
            @test go_view isa OligView{GappedOlig{DegenerateOlig}}
            @test String(go_view) == "C--GN-"
            @test n_deg_pos(go_view) == 1
            @test n_unique_oligs(go_view) == 4
            @test collect(nondegens(go_view)) == [Olig("C--GA-"), Olig("C--GC-"), Olig("C--GG-"), Olig("C--GT-")]
            
            # Gaps at start/end
            go_start = GappedOlig(olig, [1=>2])
            @test String(go_start) == "--ACGTACGT"
            @test length(go_start) == 10
            @test go_start[1] == '-' && go_start[3] == 'A'
            go_end = GappedOlig(olig, [9=>2])
            @test String(go_end) == "ACGTACGT--"
            @test length(go_end) == 10
            @test go_end[8] == 'T' && go_end[9] == '-'
            
            # Complex slicing
            nested_view = go[2:8][2:5]
            @test nested_view isa OligView{GappedOlig{Olig}}
            @test String(nested_view) == "C-GT"
            @test String(go[4:5]) == "-G"
            @test String(go[6:8]) == "T-A"
            
            # Concatenation edge cases
            @test_throws ErrorException go * go
            go_view = go[1:3]
            @test String(go_view * Olig("CG")) == "AC-CG"
            @test (go_view * Olig("CG")) isa Olig
            go_deg_view = go_deg[1:3]
            @test String(go_deg_view * Olig("CG")) == "AC-CG"
            @test (go_deg_view * Olig("CG")) isa Olig
            
            # SeqFold methods
            @test String(SeqFold.revcomp(go)) == "ACGT-AC--GT"
            @test SeqFold.revcomp(go) isa GappedOlig{Olig}
            @test SeqFold.revcomp(go).gaps == [3=>1, 5=>2]
            @test String(SeqFold.complement(go)) == "TG--CA-TGCA"
            @test SeqFold.complement(go).gaps == [3=>2, 5=>1]
            @test SeqFold.gc_content(go) ≈ SeqFold.gc_content(olig)
            @test_throws ErrorException SeqFold.fold(go)
            @test_throws ErrorException SeqFold.dg(go)
            @test_throws ErrorException SeqFold.tm(go)
            
            # Iteration
            @test collect(go) == ['A', 'C', '-', '-', 'G', 'T', '-', 'A', 'C', 'G', 'T']
            @test collect(go_deg) == ['A', 'C', '-', '-', 'G', 'N', '-', 'R', 'T']
            empty_iter = GappedOlig(Olig(""), [])
            @test iterate(empty_iter) === nothing
            
            # Nondegens for degenerate parent
            @test collect(nondegens(go_deg)) == [go_deg]
            @test collect(nondegens(go)) == [go]
        end

        # @testset "Extra" begin
        #     @testset "GappedOlig with DegenerateOlig" begin
        #         deg_olig = DegenerateOlig("ACGNRT", "deg gap test")
        #         gaps = [3=>2, 6=>1]
        #         go = GappedOlig(deg_olig, gaps)
        #         @test go isa GappedOlig{DegenerateOlig}
        #         @test String(go) == "AC--GN-RT"  # 6 bases + 3 gaps = 9 length
        #         @test go.parent == deg_olig
        #         @test go.gaps == gaps
        #         @test length(go) == 9
        #         @test n_deg_pos(go) == 2  # N and R from parent
        #         @test n_unique_oligs(go) == 8  # N=4, R=2
        #         @test collect(go) == ['A','C','-','-','G','N','-','R','T']
        #         @test String(go[2:7]) == "C--GN-"  # Includes gaps and degenerate
        #         @test String(go[8:end]) == "RT"
        #         @test go[3] == '-'
        #         @test go[6] == 'N'
        #         go_view = go[3:6]
        #         @test go_view isa OligView{GappedOlig{DegenerateOlig}}
        #         @test String(go_view) == "--GN"
        #         @test n_deg_pos(go_view) == 1
        #         @test n_unique_oligs(go_view) == 4
        #         @test collect(nondegens(go_view)) == [Olig("--GA"), Olig("--GC"), Olig("--GG"), Olig("--GT")]

        #         deg_olig = DegenerateOlig("ACGNRT", "deg gap test")
        #         # Gaps at start
        #         go_start = GappedOlig(deg_olig, [1=>2])
        #         @test String(go_start) == "--ACGNRT"
        #         @test length(go_start) == 8
        #         @test go_start[1] == '-' && go_start[2] == '-' && go_start[3] == 'A'
        #         # Gaps at end
        #         go_end = GappedOlig(deg_olig, [7=>2])
        #         @test String(go_end) == "ACGNRT--"
        #         @test length(go_end) == 8
        #         @test go_end[6] == 'T' && go_end[7] == '-' && go_end[8] == '-'
        #         # Multiple gaps with degenerate bases
        #         go_multi = GappedOlig(deg_olig, [2=>1, 4=>1, 6=>1])
        #         @test String(go_multi) == "A-C-G-N-RT"
        #         @test length(go_multi) == 9
        #         @test collect(go_multi) == ['A', '-', 'C', '-', 'G', '-', 'N', 'R', 'T']
        #         # Slice including partial gap
        #         @test String(go_multi[2:5]) == "-C-G"
        #         # Gap-only slice
        #         go = GappedOlig(deg_olig, [3=>2, 6=>1])
        #         @test String(go[3:4]) == "--"
        #     end

        #     @testset "Random Sampling" begin
        #         deg_olig = DegenerateOlig("ACRN", "rand test")
        #         rng = Random.MersenneTwister(42)
        #         rand_olig = rand(rng, deg_olig)
        #         @test rand_olig isa Olig
        #         @test String(rand_olig) in ["ACAA", "ACAG", "ACGA", "ACGG"]
        #         @test description(rand_olig) == "rand test"
        #         deg_view = deg_olig[2:4]
        #         rand_view = rand(rng, deg_view)
        #         @test rand_view isa Olig
        #         @test String(rand_view) in ["CAA", "CAG", "CGA", "CGG"]
        #         @test description(rand_view) == "rand test"
        #         # Test multiple samples to ensure coverage
        #         samples = Set(String(rand(rng, deg_olig)) for _ in 1:100)
        #         @test length(samples) > 1  # Likely to hit multiple variants
        #     end

        #     @testset "Complex Slicing for GappedOlig" begin
        #         olig = Olig("ACGTACGT", "slice test")
        #         go = GappedOlig(olig, [3=>2, 7=>1])
        #         @test String(go) == "AC--GT-ACGT"
        #         @test String(go[3:4]) == "--"  # Only gaps
        #         @test String(go[4:5]) == "-G"  # Partial gap
        #         @test String(go[6:8]) == "T-A"  # Gap in middle
        #         @test String(go[1:2]) == "AC"  # Before gaps
        #         @test String(go[end:end]) == "T"  # Single position
        #         nested_view = go[2:8][2:5]
        #         @test nested_view isa OligView{GappedOlig{Olig}}
        #         @test String(nested_view) == "--GT"  # Nested slice
        #         deg_olig = DegenerateOlig("ACGNRT", "deg slice")
        #         go_deg = GappedOlig(deg_olig, [3=>2])
        #         @test String(go_deg[2:5]) == "C--G"
        #         @test String(go_deg[3:4]) == "--"
        #     end

        #     @testset "Concatenation Edge Cases" begin
        #         go = GappedOlig(Olig("ACGT", "test"), [3=>1])
        #         @test_throws ErrorException go * go  # Gapped concatenation not supported
        #         go_view = go[1:3]
        #         @test go_view isa OligView{GappedOlig{Olig}}
        #         @test String(go_view * Olig("CG")) == "AC-CG"
        #         @test (go_view * Olig("CG")) isa Olig
        #         deg_olig = DegenerateOlig("ACN")
        #         go_deg = GappedOlig(deg_olig, [2=>1])
        #         go_deg_view = go_deg[1:3]
        #         @test String(go_deg_view * Olig("CG")) == "A-CN-CG"
        #         @test (go_deg_view * Olig("CG")) isa DegenerateOlig
        #         @test n_deg_pos(go_deg_view * Olig("CG")) == 1
        #         @test n_unique_oligs(go_deg_view * Olig("CG")) == 4
        #     end

        #     @testset "SeqFold Methods" begin
        #         deg_olig = DegenerateOlig("ACGN", "tm test")
        #         tm_result = SeqFold.tm(deg_olig, conditions=:pcr)
        #         @test tm_result.mean isa Float64
        #         @test tm_result.conf[1] <= tm_result.mean <= tm_result.conf[2]
        #         @test all(x in ["ACGA", "ACGC", "ACGG", "ACGT"] for x in nondegens(deg_olig))
        #         deg_view = deg_olig[2:4]
        #         tm_view = SeqFold.tm(deg_view, conditions=:pcr)
        #         @test tm_view.mean isa Float64
        #         @test tm_view.conf[1] <= tm_view.mean <= tm_view.conf[2]
        #         @test SeqFold.gc_content(deg_olig) ≈ (0.0 + 1.0 + 0.5 + 0.5) / 4  # A,C,G,N weights
        #         @test SeqFold.gc_content(deg_view) ≈ (1.0 + 0.5 + 0.0) / 3  # C,G,N
        #         no_gap_go = GappedOlig(Olig("ACGT"), Pair{Int,Int}[])
        #         @test SeqFold.tm(no_gap_go, conditions=:pcr).mean == SeqFold.tm("ACGT", conditions=:pcr)
        #     end
            
        # end
    end
end
