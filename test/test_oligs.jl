import AlignedOligs.Oligs: NON_DEGEN_BASES, DEGEN_BASES, ALL_BASES, IUPAC_COUNTS

const NUM_RANDOM_TESTS = 10

rseq(len, bases) = join(rand(bases, len))
rdesc() = join(rand('1':'z', rand(0:150)))
rolg(::Type{Olig}, len) =  Olig(rseq(len, NON_DEGEN_BASES), rdesc())
rolg(::Type{DegenOlig}, len) =  DegenOlig(rseq(len, ALL_BASES), rdesc())
function rolg(::Type{GappedOlig}, len)
    len == 0 && return GappedOlig("", rdesc())
    seq_chars = collect(ALL_BASES)
    append!(seq_chars, fill('-', 5)) 
    seq = rseq(len, seq_chars)
    while count(==('-'), seq) == len
        seq = rseq(len, seq_chars)
    end
    return GappedOlig(seq, rdesc())
end
@testset "Types" begin
    @test Olig <: AbstractOlig <: AbstractString
    @test DegenOlig <: AbstractDegen <: AbstractOlig
    @test GappedOlig <: AbstractGapped <: AbstractDegen <: AbstractOlig

    @test OligView <: AbstractOlig
    @test OligView{Union{Olig, DegenOlig, GappedOlig}} <: AbstractOlig
    @test_throws TypeError OligView{Int}
    
    @test NonDegenIterator{Union{Olig, DegenOlig, GappedOlig}} <: NonDegenIterator
    @test_throws TypeError NonDegenIterator{Int}
end

@testset "Olig Construction" begin
    @test Olig() == Olig("", "") == Olig("")
    @test Olig() !== Olig("", "123")
    
    for _ in 1:NUM_RANDOM_TESTS
        len = rand(0:50)
        seq = rseq(len, NON_DEGEN_BASES)
        descr = rdesc()
        olig = Olig(seq, descr)
        @test olig isa Olig
        @test String(olig) == uppercase(seq)
        @test description(olig) == descr
        
        @test Olig(DegenOlig(seq)) == olig
        @test Olig(GappedOlig(seq)) == olig
    end
    
    for bad_char in setdiff('A':'Z', NON_DEGEN_BASES)
        @test_throws ErrorException Olig("ACGT$bad_char")
    end
end

@testset "DegenOlig Construction" begin
    @test DegenOlig() == DegenOlig("")
    
    for _ in 1:NUM_RANDOM_TESTS
        len = rand(0:50)
        seq = rseq(len, ALL_BASES)
        descr = rdesc()
        deg_olig = DegenOlig(seq, descr)
        @test deg_olig isa DegenOlig
        @test String(deg_olig) == uppercase(seq)
        @test description(deg_olig) == descr
        
        if !any(c -> c in DEGEN_BASES, seq)
            nondeg = Olig(seq, descr)
            @test DegenOlig(nondeg) == deg_olig
        end
    end
    
    for bad_char in setdiff('A':'Z', ALL_BASES)
        @test_throws ErrorException DegenOlig("ACGTN$bad_char")
    end
end

@testset "GappedOlig Construction" begin
    @test GappedOlig() == GappedOlig("")
    
    for seq in ["-", "--A--", "A--", "--A", "A---T", "-G-C-", "ACGT", ""]
        go = GappedOlig(seq)
        @test String(go) == seq
        @test String(parent(go)) == filter(!=('-'), seq)
    end
    
    @test GappedOlig("A---T").gaps == [2 => 3]
    @test GappedOlig("-G-C-").gaps == [1 => 1, 2 => 1, 3 => 1]
    @test GappedOlig("ACGT").gaps == Pair{Int,Int}[]
    
    for _ in 1:NUM_RANDOM_TESTS
        len = rand(1:50)
        go = rolg(GappedOlig, len)
        seq_str = String(go)
        
        @test length(go) == len
        @test String(go) == seq_str
        @test String(parent(go)) == filter(!=('-'), seq_str)
        @test hasgaps(go) == any(==('-'), seq_str)
    end
end

@testset "Indexing and Views" begin
    for T in (Olig, DegenOlig, GappedOlig)
        for _ in 1:NUM_RANDOM_TESTS
            len = rand(0:50)
            olig = rolg(T, len)
            
            @test_throws BoundsError olig[0]

            if len > 0
                idx = rand(1:len)
                @test olig[idx] == String(olig)[idx]
                @test olig[idx] isa Char

                full_view = olig[1:len]
                @test full_view isa OligView{T}
                @test String(full_view) == String(olig)
                @test length(full_view) == len
            end
            
            if len >= 2
                start = rand(1:len-1)
                stop = rand(start:len)
                view = olig[start:stop]
                @test view isa OligView
                @test String(view) == String(olig)[start:stop]
                @test length(view) == stop - start + 1
                @test olig_range(view) == start:stop
            end
            
            if len >= 4
                outer = olig[2:len-1]
                inner = outer[1:length(outer)-1]
                @test String(inner) == String(olig)[2:len-2]
            end
        end
    end
    
    gapped = GappedOlig("A-C-G-T")
    @test gapped[1] == 'A'
    @test gapped[2] == '-'
    @test gapped[3] == 'C'
    @test gapped[4] == '-'
    @test gapped[5] == 'G'
    @test gapped[6] == '-'
    @test gapped[7] == 'T'
    
    empty_olig = Olig("")
    @test_throws BoundsError empty_olig[1]
    @test_throws BoundsError empty_olig[1:1]
    non_empty = Olig("A")
    @test_throws BoundsError non_empty[2]
    @test_throws BoundsError non_empty[0]
    @test_throws BoundsError non_empty[0:1]
end

@testset "Degeneracy Properties" begin
    nondeg = Olig("ACGT")
    @test n_deg_pos(nondeg) == 0
    @test n_unique_oligs(nondeg) == 1
    @test !hasgaps(nondeg)
    
    deg = DegenOlig("ACGTNRYS")
    @test n_deg_pos(deg) == count(c -> c in DEGEN_BASES, "NRYS")
    @test n_unique_oligs(deg) == prod(IUPAC_COUNTS[c] for c in "NRYS")
    
    gapped_deg = GappedOlig("A-N-C-")
    @test n_deg_pos(gapped_deg) == 1
    @test n_unique_oligs(gapped_deg) == 4
    @test hasgaps(gapped_deg)
    
    view = deg[3:6]
    @test n_deg_pos(view) == count(c -> c in DEGEN_BASES, "GTNR")
    @test n_unique_oligs(view) == prod(IUPAC_COUNTS[c] for c in "NR")
end

@testset "Non-Degenerate Iteration" begin
    @test isempty(collect(nondegens(Olig(""))))
    
    nondeg = Olig("ACGT", "test")
    nondeg_iter = collect(nondegens(nondeg))
    @test only(nondeg_iter) == nondeg
    @test description(only(nondeg_iter)) == "test"
    
    deg = DegenOlig("ACN", "test")
    nondeg_iter = collect(nondegens(deg))
    @test length(nondeg_iter) == 4
    @test Set(String.(nondeg_iter)) == Set(["ACA", "ACC", "ACG", "ACT"])
    @test all(description(o) == "Non-degen sample from: test" for o in nondeg_iter)
    
    gapped_deg = GappedOlig("A-N-", "test")
    @test_throws ErrorException nondegens(gapped_deg)
    
    gapped_nogaps_deg = GappedOlig("AANC", "test")
    nondeg_iter = collect(nondegens(gapped_nogaps_deg))
    @test length(nondeg_iter) == 4
    expected = ["AAAC", "AAGC", "AACC", "AATC"]
    @test Set(String.(nondeg_iter)) == Set(expected)
    
    for _ in 1:NUM_RANDOM_TESTS
        len = rand(1:8)
        deg = rolg(DegenOlig, len)
        nd_iter = nondegens(deg)
        @test length(nd_iter) == n_unique_oligs(deg)
        @test length(collect(nd_iter)) == n_unique_oligs(deg)
    end
end

@testset "Sampling Functions" begin
    for T in (Olig, DegenOlig, GappedOlig)
        for _ in 1:NUM_RANDOM_TESTS
            len = rand(1:50)
            olig = rolg(T, len)
            sampled = sampleChar(olig)
            @test sampled isa Char
            @test sampled in String(olig)
        end
        @test_throws ArgumentError sampleChar(rolg(T, 0))
    end
    
    for T in (Olig, DegenOlig, GappedOlig)
        for _ in 1:NUM_RANDOM_TESTS
            len = rand(5:50)
            olig = rolg(T, len)
            view_len = rand(1:len)
            view = sampleView(olig, view_len)
            @test view isa OligView
            @test length(view) == view_len
            @test occursin(String(view), String(olig))
        end
        empty_olig = rolg(T, 0)
        @test_throws ArgumentError sampleView(empty_olig, 1)
        non_empty = rolg(T, 3)
        @test_throws ArgumentError sampleView(non_empty, 4)
    end
    
    for _ in 1:NUM_RANDOM_TESTS
        nondeg = rolg(Olig, rand(0:50))
        @test sampleNondeg(nondeg) == nondeg
        
        deg = rolg(DegenOlig, rand(1:50))
        sampled = sampleNondeg(deg)
        @test sampled isa DegenOlig
        @test length(sampled) == length(deg)
        @test all(c in NON_DEGEN_BASES for c in String(sampled))
        
        gapped_deg = rolg(GappedOlig, rand(1:50))
        sampled = sampleNondeg(gapped_deg)
        @test sampled isa GappedOlig
        @test length(sampled) == length(gapped_deg)
        @test all(c in NON_DEGEN_BASES || c == '-' for c in String(sampled))
        @test n_deg_pos(sampled) == 0
    end
end

@testset "Conversion and Promotion" begin
    @test convert(Olig, DegenOlig("ACGT")) == Olig("ACGT")
    @test convert(DegenOlig, Olig("ACGT")) == DegenOlig("ACGT")
    @test convert(GappedOlig, Olig("ACGT")) == GappedOlig("ACGT")
    
    @test promote_type(Olig, DegenOlig) == DegenOlig
    @test promote_type(Olig, GappedOlig) == GappedOlig
    @test promote_type(DegenOlig, GappedOlig) == GappedOlig
    @test promote_type(Olig, String) == Olig
    @test promote_type(DegenOlig, SubString) == DegenOlig
    
    @test Olig("A") == "A"
    @test "A" == Olig("A")
    @test DegenOlig("N") != Olig("A")
    gapped = GappedOlig("A-C")
    @test gapped == GappedOlig("A-C")
    @test gapped != GappedOlig("AC-")
end

@testset "Iteration and Base Functions" begin
    for T in (Olig, DegenOlig, GappedOlig)
        for _ in 1:NUM_RANDOM_TESTS
            len = rand(0:50)
            olig = rolg(T, len)
            iterated = collect(olig)
            @test length(iterated) == len
            @test String(olig) == join(iterated)
        end
    end
    
    gapped = GappedOlig("--A--")
    @test collect(gapped) == ['-','-','A','-','-']
    
    for T in (Olig, DegenOlig, GappedOlig)
        for _ in 1:NUM_RANDOM_TESTS
            len = rand(0:50)
            olig = rolg(T, len)
            @test isempty(olig) == (len == 0)
            @test length(olig) == len
            @test lastindex(olig) == len
            @test ncodeunits(olig) == len
        end
    end
end