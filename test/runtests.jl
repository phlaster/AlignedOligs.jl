using AlignedOligs
using Test
using Aqua
using JET
using Random
using MAFFT_jll

Random.seed!(123)


function random_olig_string(len::Int)::String
    if len < 0
        throw(ArgumentError("Length must be non-negative"))
    end
    bases = AlignedOligs.NON_DEGEN_BASES
    return join(rand(bases, len))
end

function random_degen_string(len::Int)::String
    if len < 0
        throw(ArgumentError("Length must be non-negative"))
    end
    bases = AlignedOligs.DEGEN_BASES
    return join(rand(bases, len))
end

function random_invalid_nondeg_string(len::Int)::String
    if len < 1
        throw(ArgumentError("Length must be at least 1 for invalid strings"))
    end
    valid_str = random_olig_string(len)
    pos = rand(1:len)
    invalid_chars = setdiff(collect('A':'Z'), AlignedOligs.NON_DEGEN_BASES)
    invalid_char = rand(invalid_chars)
    return valid_str[1:pos-1] * invalid_char * valid_str[pos+1:end]
end

function random_invalid_degen_string(len::Int)::String
    if len < 1
        throw(ArgumentError("Length must be at least 1 for invalid strings"))
    end
    valid_str = random_olig_string(len)
    pos = rand(1:len)
    invalid_chars = "\\-%!@#\$%^ \""
    invalid_char = rand(invalid_chars)
    return valid_str[1:pos-1] * invalid_char * valid_str[pos+1:end]
end

function random_description()::String
    return "desc_$(rand(1:1000))"
end

function random_olig_chars(len::Int)::Vector{Char}
    return collect(random_olig_string(len))
end

@testset "AlignedOligs.jl" begin
    # @testset "Code quality (Aqua.jl)" begin
    #     Aqua.test_all(AlignedOligs)
    # end
    # @testset "Code linting (JET.jl)" begin
    #     JET.test_package(AlignedOligs; target_defined_modules = true)
    # end
    
    @testset "Oligs" begin
        @testset "Olig" begin
            empty_olig = Olig()
            empty_desc = ""
            empty_seq = ""
            empty_chars = Char[]

            @testset "Construction" begin
                # Empty constructions
                @test Olig() isa Olig
                @test Olig(empty_seq) == empty_olig
                @test Olig(empty_seq, empty_desc) == empty_olig
                @test Olig(empty_chars) == empty_olig
                @test Olig(empty_chars, empty_desc) == empty_olig
                @test description(Olig(empty_chars, random_description())) == ""

                # Random valid constructions (string, chars, with/without desc)
                for _ in 1:10
                    len = rand(1:20)
                    valid_seq = random_olig_string(len)
                    valid_chars = collect(valid_seq)
                    desc = random_description()
                    num_desc = rand(1:100)

                    olig_from_seq = Olig(valid_seq)
                    olig_from_seq_desc = Olig(valid_seq, desc)
                    olig_from_chars = Olig(valid_chars)
                    olig_from_chars_desc = Olig(valid_chars, desc)
                    olig_from_seq_num_desc = Olig(valid_seq, num_desc)

                    @test olig_from_seq isa Olig
                    @test olig_from_seq_desc isa Olig
                    @test olig_from_chars isa Olig
                    @test olig_from_chars_desc isa Olig
                    @test olig_from_seq_num_desc isa Olig

                    @test String(olig_from_seq) == valid_seq
                    @test String(olig_from_seq_desc) == valid_seq
                    @test String(olig_from_chars) == valid_seq
                    @test String(olig_from_chars_desc) == valid_seq
                    @test String(olig_from_seq_num_desc) == valid_seq

                    @test description(olig_from_seq) == ""
                    @test description(olig_from_seq_desc) == desc
                    @test description(olig_from_chars) == ""
                    @test description(olig_from_chars_desc) == desc
                    @test description(olig_from_seq_num_desc) == string(num_desc)  # Converted to string

                    @test_throws MethodError Olig(valid_seq, nothing)
                    @test Olig(valid_seq, "") isa Olig
                    @test description(Olig(valid_seq, "")) == ""
                end

                # Invalid constructions (throws ErrorException)
                for _ in 1:10
                    len = rand(1:20)
                    invalid_seq = random_invalid_nondeg_string(len)
                    invalid_chars = collect(invalid_seq)
                    desc = random_description()

                    @test_throws ErrorException Olig(invalid_seq)
                    @test_throws ErrorException Olig(invalid_seq, desc)
                    @test_throws ErrorException Olig(invalid_chars)
                    @test_throws ErrorException Olig(invalid_chars, desc)
                end                
            end

            @testset "String interface" begin
                # Empty
                @test String(empty_olig) == empty_seq
                @test length(empty_olig) == 0
                @test isempty(empty_olig)
                @test_throws BoundsError empty_olig[1]
                @test collect(empty_olig) == empty_chars

                # Random valid
                for _ in 1:10
                    len = rand(1:20)
                    valid_seq = random_olig_string(len)
                    desc = random_description()
                    olig = Olig(valid_seq, desc)

                    @test String(olig) == valid_seq
                    @test length(olig) == len
                    @test !isempty(olig)
                    @test collect(olig) == collect(valid_seq)

                    # Indexing
                    if len > 0
                        @test olig[1] == valid_seq[1]
                        @test olig[end] == valid_seq[end]
                        mid_start = rand(1:len-1)
                        mid_end = rand(mid_start+1:len)
                        @test olig[mid_start:mid_end] == valid_seq[mid_start:mid_end]
                        @test olig[end:-1:1] == reverse(valid_seq)

                        @test first(olig) == valid_seq[1]
                        @test last(olig) == valid_seq[end]

                        # occursin
                        sub_len = rand(1:len)
                        sub_start = rand(1:len-sub_len+1)
                        sub = valid_seq[sub_start:sub_start+sub_len-1]
                        @test occursin(sub, olig)
                        invalid_sub = random_invalid_nondeg_string(sub_len)
                        @test !occursin(invalid_sub, olig)
                    end

                    # Bounds errors
                    @test_throws BoundsError olig[0]
                    @test_throws BoundsError olig[len+1]
                    @test_throws BoundsError olig[1:len+1]
                end
            end

            @testset "Description" begin
                # Empty
                @test description(empty_olig) == ""

                # Random
                for _ in 1:10
                    len = rand(1:20)
                    valid_seq = random_olig_string(len)
                    desc = random_description()
                    olig_no_desc = Olig(valid_seq)
                    olig_with_desc = Olig(valid_seq, desc)

                    @test description(olig_no_desc) == ""
                    @test description(olig_with_desc) == desc
                    @test description(Olig(valid_seq, "")) == ""
                end
            end

            @testset "Conversion" begin
                # From DegenerateOlig (only non-degen)
                for _ in 1:10
                    len = rand(1:20)
                    valid_seq = random_olig_string(len)
                    desc = random_description()
                    deg_olig = DegenerateOlig(valid_seq, desc)
                    converted_olig = Olig(deg_olig)

                    @test converted_olig isa Olig
                    @test String(converted_olig) == valid_seq
                    @test description(converted_olig) == desc

                    @test convert(Olig, deg_olig) isa Olig
                    @test convert(Olig, deg_olig) == converted_olig
                end

                # Throws on degenerate
                degen_seq = "ACGN"
                deg_olig_with_degen = DegenerateOlig(degen_seq)
                @test_throws InexactError Olig(deg_olig_with_degen)
                @test_throws InexactError convert(Olig, deg_olig_with_degen)
            end

            @testset "Concatenation" begin
                # Empty
                @test empty_olig * empty_olig == empty_olig

                # Random
                for _ in 1:10
                    len1 = rand(1:10)
                    len2 = rand(1:10)
                    seq1 = random_olig_string(len1)
                    seq2 = random_olig_string(len2)
                    desc1 = random_description()
                    desc2 = random_description()
                    olig1 = Olig(seq1, desc1)
                    olig2 = Olig(seq2, desc2)

                    concat_olig = olig1 * olig2
                    @test concat_olig isa Olig
                    @test String(concat_olig) == seq1 * seq2
                    @test description(concat_olig) == "concat"
                    @test length(concat_olig) == len1 + len2

                    # Multi-concat
                    olig_a = Olig("A")
                    olig_c = Olig("C")
                    olig_g = Olig("G")
                    olig_t = Olig("T")
                    multi_concat = olig_a * olig_c * olig_g * olig_t
                    @test String(multi_concat) == "ACGT"
                end
            end

            @testset "Equality" begin
                # Empty
                @test empty_olig == empty_olig
                @test empty_olig == empty_seq
                @test empty_seq == empty_olig

                # Random
                for _ in 1:10
                    len = rand(1:20)
                    seq1 = random_olig_string(len)
                    seq2 = random_olig_string(len)
                    while seq1 == seq2  # Ensure different
                        seq2 = random_olig_string(len)
                    end
                    olig1 = Olig(seq1)
                    olig2 = Olig(seq2)
                    olig1_copy = Olig(seq1)

                    @test olig1 == olig1_copy
                    @test olig1 != olig2
                    @test olig1 == seq1
                    @test seq1 == olig1
                    @test olig1 != seq2
                    @test seq2 != olig1
                end

                # Description doesn't affect equality
                desc1 = random_description()
                desc2 = random_description()
                same_seq = random_olig_string(5)
                olig_diff_desc1 = Olig(same_seq, desc1)
                olig_diff_desc2 = Olig(same_seq, desc2)
                @test olig_diff_desc1 == olig_diff_desc2
            end

            @testset "Case handling" begin
                for _ in 1:10
                    len = rand(1:20)
                    lower_seq = lowercase(random_olig_string(len))
                    upper_seq = uppercase(lower_seq)
                    olig_lower = Olig(lower_seq)
                    olig_upper = Olig(upper_seq)

                    @test olig_lower == olig_upper
                    @test String(olig_lower) == upper_seq
                end
            end

            @testset "Empty sequence" begin
                empty_olig_with_desc = Olig(empty_seq, "empty")

                @test isempty(empty_olig_with_desc)
                @test length(empty_olig_with_desc) == 0
                @test String(empty_olig_with_desc) == empty_seq
                @test_throws BoundsError empty_olig_with_desc[1]
                @test description(empty_olig_with_desc) == ""
            end

            @testset "Display" begin
                # Fixed cases for reproducibility
                short_seq = "AGTC"
                short_desc = "descr"
                no_desc = ""
                long_seq = "A"^25

                olig_short_desc = Olig(short_seq, short_desc)
                olig_short_no_desc = Olig(short_seq, no_desc)
                olig_long_no_desc = Olig(long_seq, no_desc)

                # Short display
                @test sprint(show, olig_short_desc) == "Olig(\"AGTC\", len=4, desc=\"descr\")"
                @test sprint(show, olig_short_no_desc) == "Olig(\"AGTC\", len=4)"
                @test sprint(show, olig_long_no_desc) == "Olig(\"AAAAAAAAAAAAAAAAA...\", len=25)"

                # Full display
                @test sprint(show, MIME"text/plain"(), olig_short_desc) == "Olig\n  Sequence: AGTC\n  Length: 4\n  Description: \"descr\""
                @test sprint(show, MIME"text/plain"(), olig_short_no_desc) == "Olig\n  Sequence: AGTC\n  Length: 4\n  Description: (none)"
                @test sprint(show, MIME"text/plain"(), olig_long_no_desc) == "Olig\n  Sequence: AAAAAAAAAAAAAAAAAAAAAAAAA\n  Length: 25\n  Description: (none)"

                # Random short and long
                for _ in 1:5
                    rand_len_short = rand(1:20)
                    rand_seq_short = random_olig_string(rand_len_short)
                    rand_desc_short = random_description()
                    olig_rand_short = Olig(rand_seq_short, rand_desc_short)

                    expected_short_show = "Olig(\"$rand_seq_short\", len=$rand_len_short, desc=\"$rand_desc_short\")"
                    @test sprint(show, olig_rand_short) == expected_short_show

                    expected_full_show = "Olig\n  Sequence: $rand_seq_short\n  Length: $rand_len_short\n  Description: \"$rand_desc_short\""
                    @test sprint(show, MIME"text/plain"(), olig_rand_short) == expected_full_show
                end

                for _ in 1:5
                    rand_len_long = rand(21:50)
                    rand_seq_long = random_olig_string(rand_len_long)
                    olig_rand_long = Olig(rand_seq_long)

                    truncated_seq = rand_seq_long[1:17] * "..."
                    expected_short_show = "Olig(\"$truncated_seq\", len=$rand_len_long)"
                    @test sprint(show, olig_rand_long) == expected_short_show

                    expected_full_show = "Olig\n  Sequence: $rand_seq_long\n  Length: $rand_len_long\n  Description: (none)"
                    @test sprint(show, MIME"text/plain"(), olig_rand_long) == expected_full_show
                end
            end
        end

        @testset "DegenerateOlig" begin
            empty_degen = DegenerateOlig("")
            empty_desc = ""
            empty_seq = ""
            empty_chars = Char[]

            @testset "Construction" begin
                # Empty constructions
                @test DegenerateOlig() == empty_degen
                @test DegenerateOlig(empty_seq) == empty_degen
                @test DegenerateOlig(empty_seq, empty_desc) == empty_degen
                @test DegenerateOlig(empty_chars) == empty_degen
                @test DegenerateOlig(empty_chars, empty_desc) == empty_degen
                @test description(DegenerateOlig(empty_chars, random_description())) == ""

                # Random valid constructions (string, chars, with/without desc)
                for _ in 1:10
                    len = rand(1:20)
                    valid_seq = random_degen_string(len)
                    valid_chars = collect(valid_seq)
                    desc = random_description()
                    num_desc = rand(1:100)

                    deg_from_seq = DegenerateOlig(valid_seq)
                    deg_from_seq_desc = DegenerateOlig(valid_seq, desc)
                    deg_from_chars = DegenerateOlig(valid_chars)
                    deg_from_chars_desc = DegenerateOlig(valid_chars, desc)
                    deg_from_seq_num_desc = DegenerateOlig(valid_seq, num_desc)

                    @test deg_from_seq isa DegenerateOlig
                    @test deg_from_seq_desc isa DegenerateOlig
                    @test deg_from_chars isa DegenerateOlig
                    @test deg_from_chars_desc isa DegenerateOlig
                    @test deg_from_seq_num_desc isa DegenerateOlig

                    @test String(deg_from_seq) == uppercase(valid_seq)
                    @test String(deg_from_seq_desc) == uppercase(valid_seq)
                    @test String(deg_from_chars) == uppercase(valid_seq)
                    @test String(deg_from_chars_desc) == uppercase(valid_seq)
                    @test String(deg_from_seq_num_desc) == uppercase(valid_seq)

                    @test description(deg_from_seq) == ""
                    @test description(deg_from_seq_desc) == desc
                    @test description(deg_from_chars) == ""
                    @test description(deg_from_chars_desc) == desc
                    @test description(deg_from_seq_num_desc) == string(num_desc)  # Converted to string

                    # Properties
                    expected_deg_pos = count(c -> c in AlignedOligs.DEGEN_BASES, uppercase(valid_seq))
                    expected_unique = reduce(*, (AlignedOligs.IUPAC_COUNTS[c] for c in uppercase(valid_seq)); init=BigInt(1))
                    @test n_deg_pos(deg_from_seq) == expected_deg_pos
                    @test n_unique_oligs(deg_from_seq) == expected_unique
                    @test n_deg_pos(deg_from_seq_desc) == expected_deg_pos
                    @test n_unique_oligs(deg_from_seq_desc) == expected_unique

                    @test_throws MethodError DegenerateOlig(valid_seq, nothing)
                    @test DegenerateOlig(valid_seq, "") isa DegenerateOlig
                    @test description(DegenerateOlig(valid_seq, "")) == ""
                end

                # Invalid constructions (throws ErrorException)
                for _ in 1:10
                    len = rand(1:20)
                    invalid_seq = random_invalid_degen_string(len)
                    invalid_chars = collect(invalid_seq)
                    desc = random_description()

                    @test_throws ErrorException DegenerateOlig(invalid_seq)
                    @test_throws ErrorException DegenerateOlig(invalid_seq, desc)
                    @test_throws ErrorException DegenerateOlig(invalid_chars)
                    @test_throws ErrorException DegenerateOlig(invalid_chars, desc)
                end                
            end

            @testset "String interface" begin
                # Empty
                @test String(empty_degen) == empty_seq
                @test length(empty_degen) == 0
                @test isempty(empty_degen)
                @test_throws BoundsError empty_degen[1]
                @test collect(empty_degen) == empty_chars

                # Random valid
                for _ in 1:10
                    len = rand(1:20)
                    valid_seq = random_degen_string(len)
                    desc = random_description()
                    deg = DegenerateOlig(valid_seq, desc)

                    @test String(deg) == uppercase(valid_seq)
                    @test length(deg) == len
                    @test !isempty(deg)
                    @test collect(deg) == collect(uppercase(valid_seq))

                    # Indexing
                    if len > 0
                        @test deg[1] == uppercase(valid_seq)[1]
                        @test deg[end] == uppercase(valid_seq)[end]
                        mid_start = rand(1:len-1)
                        mid_end = rand(mid_start+1:len)
                        @test deg[mid_start:mid_end] == uppercase(valid_seq)[mid_start:mid_end]
                        @test deg[end:-1:1] == reverse(uppercase(valid_seq))

                        @test first(deg) == uppercase(valid_seq)[1]
                        @test last(deg) == uppercase(valid_seq)[end]

                        # occursin
                        sub_len = rand(1:len)
                        sub_start = rand(1:len-sub_len+1)
                        sub = uppercase(valid_seq)[sub_start:sub_start+sub_len-1]
                        @test occursin(sub, deg)
                        invalid_sub = random_invalid_degen_string(sub_len)
                        @test !occursin(invalid_sub, deg)
                    end

                    # Bounds errors
                    @test_throws BoundsError deg[0]
                    @test_throws BoundsError deg[len+1]
                    @test_throws BoundsError deg[1:len+1]
                end
            end

            @testset "Description" begin
                # Empty
                @test description(empty_degen) == ""

                # Random
                for _ in 1:10
                    len = rand(1:20)
                    valid_seq = random_degen_string(len)
                    desc = random_description()
                    deg_no_desc = DegenerateOlig(valid_seq)
                    deg_with_desc = DegenerateOlig(valid_seq, desc)

                    @test description(deg_no_desc) == ""
                    @test description(deg_with_desc) == desc
                    @test description(DegenerateOlig(valid_seq, "")) == ""
                end
            end

            @testset "Properties" begin
                # Specific cases
                deg_olig = DegenerateOlig("ACGN", "test")
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

                # Empty
                @test n_deg_pos(empty_degen) == 0
                @test n_unique_oligs(empty_degen) == 1

                # Random
                for _ in 1:10
                    len = rand(1:20)
                    valid_seq = random_degen_string(len)
                    deg = DegenerateOlig(valid_seq)
                    expected_deg_pos = count(c -> c in random_degen_bases(), uppercase(valid_seq))
                    expected_unique = reduce(*, (AlignedOligs.IUPAC_COUNTS[c] for c in uppercase(valid_seq)); init=BigInt(1))
                    @test n_deg_pos(deg) == expected_deg_pos
                    @test n_unique_oligs(deg) == expected_unique
                end
            end

            @testset "Conversion" begin
                # From Olig
                for _ in 1:10
                    len = rand(1:20)
                    valid_seq = random_olig_string(len)  # non-degen
                    desc = random_description()
                    olig = Olig(valid_seq, desc)
                    converted_deg = DegenerateOlig(olig)

                    @test converted_deg isa DegenerateOlig
                    @test String(converted_deg) == valid_seq
                    @test description(converted_deg) == desc
                    @test n_deg_pos(converted_deg) == 0
                    @test n_unique_oligs(converted_deg) == 1

                    @test convert(DegenerateOlig, olig) isa DegenerateOlig
                    @test convert(DegenerateOlig, olig) == converted_deg
                end

                # From DegenerateOlig with degen (self)
                degen_seq = "ACGN"
                deg_olig_with_degen = DegenerateOlig(degen_seq)
                converted_self = DegenerateOlig(deg_olig_with_degen)
                @test converted_self == deg_olig_with_degen
                @test n_deg_pos(converted_self) == 1
                @test n_unique_oligs(converted_self) == 4
            end

            @testset "Concatenation" begin
                # Empty
                @test empty_degen * empty_degen == empty_degen

                # Random
                for _ in 1:10
                    len1 = rand(1:10)
                    len2 = rand(1:10)
                    seq1 = random_degen_string(len1)
                    seq2 = random_degen_string(len2)
                    desc1 = random_description()
                    desc2 = random_description()
                    deg1 = DegenerateOlig(seq1, desc1)
                    deg2 = DegenerateOlig(seq2, desc2)

                    concat_deg = deg1 * deg2
                    expected_type = (n_deg_pos(deg1) > 0 || n_deg_pos(deg2) > 0) ? DegenerateOlig : Olig
                    @test concat_deg isa expected_type
                    @test String(concat_deg) == uppercase(seq1) * uppercase(seq2)
                    @test description(concat_deg) == "concat"
                    @test length(concat_deg) == len1 + len2

                    expected_deg_pos = n_deg_pos(deg1) + n_deg_pos(deg2)
                    expected_unique = n_unique_oligs(deg1) * n_unique_oligs(deg2)
                    if concat_deg isa DegenerateOlig
                        @test n_deg_pos(concat_deg) == expected_deg_pos
                        @test n_unique_oligs(concat_deg) == expected_unique
                    end
                end

                # Multi-concat with mixed
                deg_a = DegenerateOlig("A")
                deg_c = DegenerateOlig("CN")
                deg_g = DegenerateOlig("G")
                deg_t = DegenerateOlig("T")
                multi_concat = deg_a * deg_c * deg_g * deg_t
                @test multi_concat isa DegenerateOlig
                @test String(multi_concat) == "ACNGT"
            end

            @testset "Equality" begin
                # Empty
                @test empty_degen == empty_degen
                @test empty_degen == empty_seq
                @test empty_seq == empty_degen

                # Random
                for _ in 1:10
                    len = rand(1:20)
                    seq1 = random_degen_string(len)
                    seq2 = random_degen_string(len)
                    while uppercase(seq1) == uppercase(seq2)  # Ensure different
                        seq2 = random_degen_string(len)
                    end
                    deg1 = DegenerateOlig(seq1)
                    deg2 = DegenerateOlig(seq2)
                    deg1_copy = DegenerateOlig(seq1)

                    @test deg1 == deg1_copy
                    @test deg1 != deg2
                    @test deg1 == uppercase(seq1)
                    @test uppercase(seq1) == deg1
                    @test deg1 != uppercase(seq2)
                    @test uppercase(seq2) != deg1
                end

                # Description doesn't affect equality
                desc1 = random_description()
                desc2 = random_description()
                same_seq = random_degen_string(5)
                deg_diff_desc1 = DegenerateOlig(same_seq, desc1)
                deg_diff_desc2 = DegenerateOlig(same_seq, desc2)
                @test deg_diff_desc1 == deg_diff_desc2
            end

            @testset "Case handling" begin
                for _ in 1:10
                    len = rand(1:20)
                    lower_seq = lowercase(random_degen_string(len))
                    upper_seq = uppercase(lower_seq)
                    deg_lower = DegenerateOlig(lower_seq)
                    deg_upper = DegenerateOlig(upper_seq)

                    @test deg_lower == deg_upper
                    @test String(deg_lower) == upper_seq
                end
            end

            @testset "Empty sequence" begin
                empty_degen_with_desc = DegenerateOlig(empty_seq, "empty")

                @test isempty(empty_degen_with_desc)
                @test length(empty_degen_with_desc) == 0
                @test String(empty_degen_with_desc) == empty_seq
                @test_throws BoundsError empty_degen_with_desc[1]
                @test description(empty_degen_with_desc) == ""
            end

            @testset "Display" begin
                # Fixed cases for reproducibility
                short_seq = "ACGN"
                short_desc = "descr"
                no_desc = ""
                long_seq = "N"^25

                deg_short_desc = DegenerateOlig(short_seq, short_desc)
                deg_short_no_desc = DegenerateOlig(short_seq, no_desc)
                deg_long_no_desc = DegenerateOlig(long_seq, no_desc)

                # Short display
                @test sprint(show, deg_short_desc) == "DegenerateOlig(\"ACGN\", len=4, desc=\"descr\")"
                @test sprint(show, deg_short_no_desc) == "DegenerateOlig(\"ACGN\", len=4)"
                @test sprint(show, deg_long_no_desc) == "DegenerateOlig(\"NNNNNNNNNNNNNNNNN...\", len=25)"

                # Full display
                @test sprint(show, MIME"text/plain"(), deg_short_desc) == "DegenerateOlig\n  Sequence: ACGN\n  Length: 4\n  Description: \"descr\""
                @test sprint(show, MIME"text/plain"(), deg_short_no_desc) == "DegenerateOlig\n  Sequence: ACGN\n  Length: 4\n  Description: (none)"
                @test sprint(show, MIME"text/plain"(), deg_long_no_desc) == "DegenerateOlig\n  Sequence: NNNNNNNNNNNNNNNNNNNNNNNNN\n  Length: 25\n  Description: (none)"

                # Random short and long
                for _ in 1:5
                    rand_len_short = rand(1:20)
                    rand_seq_short = random_degen_string(rand_len_short)
                    rand_desc_short = random_description()
                    deg_rand_short = DegenerateOlig(rand_seq_short, rand_desc_short)

                    expected_short_show = "DegenerateOlig(\"$(uppercase(rand_seq_short))\", len=$rand_len_short, desc=\"$rand_desc_short\")"
                    @test sprint(show, deg_rand_short) == expected_short_show

                    expected_full_show = "DegenerateOlig\n  Sequence: $(uppercase(rand_seq_short))\n  Length: $rand_len_short\n  Description: \"$rand_desc_short\""
                    @test sprint(show, MIME"text/plain"(), deg_rand_short) == expected_full_show
                end

                for _ in 1:5
                    rand_len_long = rand(21:50)
                    rand_seq_long = random_degen_string(rand_len_long)
                    deg_rand_long = DegenerateOlig(rand_seq_long)

                    truncated_seq = uppercase(rand_seq_long)[1:17] * "..."
                    expected_short_show = "DegenerateOlig(\"$truncated_seq\", len=$rand_len_long)"
                    @test sprint(show, deg_rand_long) == expected_short_show

                    expected_full_show = "DegenerateOlig\n  Sequence: $(uppercase(rand_seq_long))\n  Length: $rand_len_long\n  Description: (none)"
                    @test sprint(show, MIME"text/plain"(), deg_rand_long) == expected_full_show
                end
            end
        end
    end
end