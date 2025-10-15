const PRIMER_LEN_RANGE=18:22
const GC_RANGE=40:60
const TM_RANGE=55:60
const NONDEGEN_PRIMER_TAIL=3
const HEAD_DEGEN=5
const HEAD_SLACK=0.15
const TAIL_DEGEN=0
const TAIL_SLACK=0.05
const MINIMAL_DELTA_G=-5.0 # kcal/mol
const MAXIMAL_GAP_PROPORTION=0.25
const MAX_NONDEG_OLIG_VARIANTS=100

import FastaIO

abstract type AbstractOlig <: AbstractString end

Base.ncodeunits(olig::AbstractOlig) = ncodeunits(olig.seq)
Base.codeunit(olig::AbstractOlig) = codeunit(olig.seq)
Base.codeunit(olig::AbstractOlig, i::Integer) = codeunit(olig.seq, i)
Base.isvalid(olig::AbstractOlig, i::Integer) = isvalid(olig.seq, i)
Base.iterate(olig::AbstractOlig) = iterate(olig.seq)
Base.iterate(olig::AbstractOlig, state::Integer) = iterate(olig.seq, state)
Base.lastindex(olig::AbstractOlig) = lastindex(olig.seq)
Base.length(olig::AbstractOlig) = length(olig.seq)
Base.isempty(olig::AbstractOlig) = isempty(olig.seq)

struct Olig <: AbstractOlig
    seq::String
    
    function Olig(seq::AbstractString)
        if isempty(seq)
            return new("")
        end
        seq = uppercase(seq)
        allowed_chars = Set("ACGT")
        seq_chars = Set(seq)
        @assert allowed_chars >= seq_chars "Olig contains unallowed characters: $(join(setdiff(seq_chars, allowed_chars) |> collect, ", "))"
        return new(seq)
    end
    
end

const EMPTY_OLIG = Olig("")

Olig() = EMPTY_OLIG
Olig(chars::Vector{Char}) = Olig(String(chars))

function Base.show(io::IO, olig::Olig)
    max_width = 20
    if length(olig.seq) > max_width
        seq_display = olig.seq[1:min(end, max_width-3)] * "..."
    else
        seq_display = olig.seq
    end
    print(io, "Olig(\"", seq_display, "\", len=$(length(olig)))")
end

n_degens(d::Olig) = 0
n_unique_oligs(d::Olig) = 1
revcomp(olig::Olig) = Olig(revcomp(olig.seq))
complement(olig::Olig) = Olig(complement(olig.seq))

struct DegenerateOlig <: AbstractOlig
    seq::String
    n_deg_pos::Int
    n_unique_oligs::BigInt


    function DegenerateOlig(seq::AbstractString)
        if isempty(seq)
            return new("", 0, 1)
        end

        seq = uppercase(seq)
        allowed_chars = Set("ACGTMRWSYKVHDBN")
        seq_chars = Set(seq)
        @assert allowed_chars >= seq_chars "DegenerateOlig contains unallowed characters: $(join(setdiff(seq_chars, allowed_chars) |> collect, ", "))"
        
        n_degenerate = 0
        n_possible = BigInt(1)
        
        iupac_counts = Dict(
            'A' => 1, 'C' => 1, 'G' => 1, 'T' => 1,
            'M' => 2, 'R' => 2, 'W' => 2, 'S' => 2, 'Y' => 2, 'K' => 2,
            'V' => 3, 'H' => 3, 'D' => 3, 'B' => 3,
            'N' => 4
        )
        
        for char in seq
            if char in "MRWSYKVHDBN"
                n_degenerate += 1
            end
            n_possible *= iupac_counts[char]
        end 
        return new(seq, n_degenerate, n_possible)
    end
    

    Base.:*(deg::DegenerateOlig, olig::Olig) = new(deg.seq * olig.seq, deg.n_deg_pos, deg.n_unique_oligs)
    Base.:*(olig::Olig, deg::DegenerateOlig) = new(olig.seq * deg.seq, deg.n_deg_pos, deg.n_unique_oligs)
    function Base.:*(d1::DegenerateOlig, d2::DegenerateOlig)
        seq = d1.seq * d2.seq
        degen_pos = d1.n_deg_pos + d2.n_deg_pos
        n_possible = d1.n_unique_oligs * d2.n_unique_oligs
        new(seq, degen_pos, n_possible)
    end
    DegenerateOlig(olig::Olig) = new(olig.seq, 0, 1)

end

const EMPTY_DEGENERATE_OLIG = DegenerateOlig("")

DegenerateOlig() = EMPTY_DEGENERATE_OLIG

DegenerateOlig(chars::Vector{Char}) = DegenerateOlig(String(chars))
Base.:*(o1::Olig, o2::Olig) = Olig(o1.seq * o2.seq)

Base.:(==)(o1::AbstractOlig, o2::AbstractOlig) = o1.seq == o2.seq

function Base.show(io::IO, olig::DegenerateOlig)
    max_width = 15
    if length(olig.seq) > max_width
        seq_display = SubString(olig.seq, 1:min(lastindex(olig.seq), max_width-3)) * "..."
    else
        seq_display = olig.seq
    end
    print(io, "DegenerateOlig(\"", seq_display, "\"), len=$(length(olig)), ")
    print(io, "with $(olig.n_deg_pos) degen pos., with $(olig.n_unique_oligs > 10000 ? ">10k" : olig.n_unique_oligs) olig variants")
end

n_degens(d::DegenerateOlig) = d.n_deg_pos
n_unique_oligs(d::AbstractOlig) = d.n_unique_oligs

revcomp(deg::DegenerateOlig) = DegenerateOlig(revcomp(deg.seq))
complement(deg::DegenerateOlig) = DegenerateOlig(complement(deg.seq))

function Olig(deg_olig::DegenerateOlig)
    if deg_olig.n_deg_pos > 0
        throw(InexactError(:Olig, Olig, deg_olig))
    end
    return Olig(deg_olig.seq)
end

function list_nondegen(olig::AbstractOlig)
    iupac_map = Dict(
        'A'=>"A",  'C'=>"C",  'G'=>"G",  'T'=>"T",
        'R'=>"AG", 'Y'=>"CT", 'S'=>"CG", 'W'=>"AT",
        'K'=>"GT", 'M'=>"AC", 'B'=>"CGT",'D'=>"AGT",
        'H'=>"ACT",'V'=>"ACG",'N'=>"ACGT"
    )
    
    primer_length = length(olig)
    options = Vector{String}(undef, primer_length)
    lens = Vector{Int}(undef, primer_length)
    
    for (i, c) in enumerate(olig)
        options[i] = iupac_map[c]
        lens[i] = length(options[i])
    end
    
    resulting_oligs = Vector{Olig}(undef, n_unique_oligs(olig))
    buffer = Vector{Char}(undef, primer_length)
    indices = ones(Int, primer_length)
    
    @inbounds for i in 1:n_unique_oligs(olig)
        @simd for j in 1:primer_length
            buffer[j] = options[j][indices[j]]
        end
        resulting_oligs[i] = Olig(copy(buffer))
        
        pos = primer_length
        while pos > 0
            indices[pos] += 1
            if indices[pos] <= lens[pos]
                break
            end
            indices[pos] = 1
            pos -= 1
        end
    end
    
    return resulting_oligs
end

list_nondegen(olig::Olig) = [olig]

Base.convert(::Type{Olig}, deg_olig::DegenerateOlig) = Olig(deg_olig)
Base.convert(::Type{DegenerateOlig}, olig::Olig) = DegenerateOlig(olig)
Base.promote_rule(::Type{Olig}, ::Type{DegenerateOlig}) = DegenerateOlig

_comp(b::UInt8)::UInt8 = b == 0x41 ? 0x54 :
    b == 0x54 ? 0x41 : b == 0x43 ? 0x47 :
    b == 0x47 ? 0x43 : b == 0x61 ? 0x74 :
    b == 0x74 ? 0x61 : b == 0x63 ? 0x67 :
    b == 0x67 ? 0x63 : b == 0x52 ? 0x59 :
    b == 0x59 ? 0x52 : b == 0x72 ? 0x79 :
    b == 0x79 ? 0x72 : b == 0x4b ? 0x4d :
    b == 0x4d ? 0x4b : b == 0x6b ? 0x6d :
    b == 0x6d ? 0x6b : b == 0x42 ? 0x56 :
    b == 0x56 ? 0x42 : b == 0x62 ? 0x76 :
    b == 0x76 ? 0x62 : b == 0x44 ? 0x48 :
    b == 0x48 ? 0x44 : b == 0x64 ? 0x68 :
    b == 0x68 ? 0x64 : b
const DNA_COMP_TABLE = collect(ntuple(i->_comp(UInt8(i-1)), Val(256)))

@inline function _revcomp_bytes!(out::Vector{UInt8}, seq::AbstractVector{UInt8})
    n = length(seq)
    @inbounds @simd for i in 1:n
        out[i] = DNA_COMP_TABLE[ seq[n - i + 1] + 1 ]
    end
    return out
end

@inline function _complement_bytes!(out::Vector{UInt8}, seq::AbstractVector{UInt8})
    n = length(seq)
    @inbounds @simd for i in 1:n
        out[i] = DNA_COMP_TABLE[ seq[i] + 1 ]
    end
    return out
end

function revcomp(seq::AbstractVector{UInt8})
    out = Vector{UInt8}(undef, length(seq))
    return _revcomp_bytes!(out, seq)
end

function revcomp(seq::AbstractString)
    out_bytes = revcomp(codeunits(seq))
    return String(out_bytes)
end

function complement(seq::AbstractVector{UInt8})
    out = Vector{UInt8}(undef, length(seq))
    return _complement_bytes!(out, seq)
end

function complement(seq::AbstractString)
    out_bytes = complement(codeunits(seq))
    return String(out_bytes)
end

function complement(nucleotide::AbstractChar)
    Char(DNA_COMP_TABLE[Int(nucleotide) + 1])
end








abstract type AbstractMSA end

struct MSA <: AbstractMSA
    sequences::Vector{String}
    base_count::Matrix{Float64}

    function _base_count_matrix(sequences::Vector{<:AbstractString})
        isempty(sequences) && return zeros(4, 0)

        aln_length = length(first(sequences))
        @assert all(==(aln_length)∘length, sequences) "All sequences in MSA must have equal length"
        basecount_mtr = zeros(4, aln_length) # A, C, G, T
        l = inv(length(sequences))
        N_scores = (0.25l, 0.25l, 0.25l, 0.25l)

        for (i, mtr_col) in enumerate(eachcol(basecount_mtr))
            for seq in sequences
                c = seq[i]
                if c == 'A'
                    mtr_col[1] += l
                elseif c == 'C'
                    mtr_col[2] += l
                elseif c == 'G'
                    mtr_col[3] += l
                elseif c == 'T'
                    mtr_col[4] += l
                elseif c == 'N' || c == '?'
                    mtr_col .+= N_scores
                elseif c == '-' || c == '.'
                    continue
                else
                    error("Wrong char found in sequence: ", seq[i])
                end
            end
        end

        return round.(basecount_mtr; digits=4)
    end

    function MSA(sequences::Vector{<:AbstractString}, gap_tolerance::Real=0.0)
        @assert 0.0 <= gap_tolerance <= 1.0 "gap_tolerance must be between 0.0 and 1.0"
        
        isempty(sequences) && return new(String[], zeros(4,0))

        _allowed_msa_letters = Set("ACGTN-?.")
        seq_length = length(first(sequences))
        @assert all(length.(sequences) .== seq_length) "Aligned sequences must have similar length"
        
        sequences = uppercase.(sequences)
        for seq in sequences
            chars_in_seq = Set(seq)
            @assert _allowed_msa_letters >= chars_in_seq "Sequence contains unallowed characters: $(join(setdiff(chars_in_seq, _allowed_msa_letters) |> collect, ", "))"
        end
        
        if length(sequences) > 0 && seq_length > 0
            filtered_sequences = fill("", length(sequences))
            
            for pos in 1:seq_length
                position_chars = [seq[pos] for seq in sequences]
                n_gaps = count(c -> c in "-?.", position_chars)
                n_sequences = length(sequences)
                n_non_gaps = n_sequences - n_gaps
                
                if n_non_gaps > 0 && (n_non_gaps / n_sequences) >= gap_tolerance
                    for seq_idx in 1:length(sequences)
                        filtered_sequences[seq_idx] *= sequences[seq_idx][pos]
                    end
                end
            end
            
            sequences = filtered_sequences
        end
        
        base_count = _base_count_matrix(sequences)
        
        return new(sequences, base_count)
    end

    function MSA(aligned_fasta::AbstractString, gap_tolerance::Real=0.0)
        data = FastaIO.readfasta(aligned_fasta)
        sequences = getindex.(data, 2)

        return MSA(sequences, gap_tolerance)
    end

    function MSA(predicate_func::Function, aligned_fasta::AbstractString, gap_tolerance::Real=0.0)
        data = FastaIO.readfasta(aligned_fasta)
        headers = getindex.(data, 1)
        valid_ids = findall(predicate_func, headers)
        sequences = getindex.(data, 2)[valid_ids]

        return MSA(sequences, gap_tolerance)
    end

    MSA() = MSA(String[])
end

struct MSAView <: AbstractMSA
    parent_msa::MSA
    interval::UnitRange{Int}
end

Base.getindex(msa::MSA, i::Int) = view(msa.base_count, :, i)
Base.getindex(msa::MSA, interval::UnitRange{Int}) = MSAView(msa, interval)

function Base.getindex(msa::MSAView, i::Int)
    checkbounds(1:length(msa), i)

    abs_index = i+msa.interval.start-1
    return view(msa.parent_msa.base_count, :, abs_index)
end

function Base.getindex(msa::MSAView, interval::UnitRange{Int})
    checkbounds(1:length(msa), interval)

    original_start = msa.interval.start
    new_start = original_start + interval.start - 1
    new_stop = original_start + interval.stop - 1
    new_range = new_start:new_stop
    return MSAView(msa.parent_msa, new_range)
end

Base.length(msa::MSA) = size(msa.base_count, 2)
Base.length(msa::MSAView) = length(msa.interval)

Base.isempty(msa::AbstractMSA) = length(msa)==0

Base.firstindex(msa::AbstractMSA) = firstindex(msa.interval)
Base.lastindex(msa::AbstractMSA) = length(msa)

nseqs(msa::MSA) = length(msa.sequences)
nseqs(msa::MSAView) = length(msa.parent_msa.sequences)

Base.size(msa::AbstractMSA) = (nseqs(msa), length(msa))
Base.size(msa::AbstractMSA, i::Int) = getindex(size(msa), i)

getrange(msa::MSA) = 1:length(msa)
getrange(msa::MSAView) = msa.interval

getsequence(msa::MSA, i::Int) = getindex(msa.sequences, i)
getsequence(msa::MSAView, i::Int) = view(msa.parent_msa.sequences[i], msa.interval)

function Base.iterate(msa::AbstractMSA, state::Int = 1)
    state > length(msa) && return nothing
    return (getindex(msa, state), state + 1)
end

get_base_count(msa::MSA) = msa.base_count
get_base_count(msa::MSAView) = view(msa.parent_msa.base_count, :, msa.interval)

dry_msa(msa::MSA) = MSA(unique(msa.sequences))
dry_msa(msa::MSAView) = MSA(unique(getindex.(msa.parent_msa.sequences, Ref(msa.interval))))

function Base.show(io::IO, msa::AbstractMSA)
    n_sequences = nseqs(msa)
    if n_sequences == 0
        print(io, typeof(msa), " with 0 sequences")
        return
    end
    seq_length = length(msa)
    print(io, typeof(msa), " with $n_sequences sequences of length $seq_length")
    if seq_length == 0
        return
    end

    local terminal_height, terminal_width
    try
        terminal_height, terminal_width = displaysize(io)
    catch
        terminal_height, terminal_width = 24, 80
    end
    max_display_height = max(5, terminal_height - 6)
    max_display_width = max(20, terminal_width - 30)

    if seq_length <= max_display_width
        println(io, ":")
        if n_sequences <= max_display_height
            for i in 1:n_sequences
                println(io, getsequence(msa, i))
            end
        else
            n_first = max_display_height ÷ 2
            n_last = max_display_height - n_first - 1
            n_first = min(n_first, n_sequences)
            n_last = min(n_last, n_sequences)
            n_last = min(n_last, n_sequences - n_first)

            for i in 1:n_first
                println(io, getsequence(msa, i))
            end
            if n_first + n_last < n_sequences
                println(io, "...")
            end
            start_idx_last = max(n_first + 1, n_sequences - n_last + 1)
            for i in start_idx_last:n_sequences
                if i > 0 && i <= n_sequences
                    println(io, getsequence(msa, i))
                end
            end
        end
    else
        println(io, ":")
        half_width = max_display_width ÷ 2 - 2
        half_width = max(1, half_width)

        if n_sequences <= max_display_height
            for i in 1:n_sequences
                seq_str = string(getsequence(msa, i))
                seq_len = length(seq_str)
                if seq_len <= max_display_width
                    println(io, seq_str)
                else
                    start_end_idx = max(1, min(half_width, seq_len))
                    end_start_idx = max(start_end_idx + 1, seq_len - half_width + 1)
                    start_part = seq_str[1:start_end_idx]
                    end_part = seq_str[end_start_idx:end]
                    println(io, start_part * "..." * end_part)
                end
            end
        else
            n_first = max_display_height ÷ 2
            n_last = max_display_height - n_first - 1
            n_first = min(n_first, n_sequences)
            n_last = min(n_last, n_sequences)
            n_last = min(n_last, n_sequences - n_first)

            for i in 1:n_first
                 seq_str = string(getsequence(msa, i))
                 seq_len = length(seq_str)
                 if seq_len <= max_display_width
                     println(io, seq_str)
                 else
                     start_end_idx = max(1, min(half_width, seq_len))
                     end_start_idx = max(start_end_idx + 1, seq_len - half_width + 1)
                     start_part = seq_str[1:start_end_idx]
                     end_part = seq_str[end_start_idx:end]
                     println(io, start_part * "..." * end_part)
                 end
            end

            if n_first + n_last < n_sequences
                println(io, " " ^ half_width * " ⋮ ")
            end

            start_idx_last = max(n_first + 1, n_sequences - n_last + 1)
            for i in start_idx_last:n_sequences
                seq_str = string(getsequence(msa, i))
                seq_len = length(seq_str)
                if seq_len <= max_display_width
                    println(io, seq_str)
                else
                    start_end_idx = max(1, min(half_width, seq_len))
                    end_start_idx = max(start_end_idx + 1, seq_len - half_width + 1)
                    start_part = seq_str[1:start_end_idx]
                    end_part = seq_str[end_start_idx:end]
                    println(io, start_part * "..." * end_part)
                end
            end
        end
    end
end

function consensus_major(msa::AbstractMSA, slack::Real=0.0)
    @assert zero(slack) <= slack <= one(slack) "slack must be between 0.0 and 1.0"
    n_positions = length(msa)
    if n_positions == 0
        return DegenerateOlig()
    end
    
    consensus_chars = Vector{Char}(undef, n_positions)
    bases = ('A', 'C', 'G', 'T')
    
    for (i, col) in enumerate(msa)
        max_freq = maximum(col)
        
        if max_freq < slack
            consensus_chars[i] = 'N'
        else
            max_idx = argmax(col)
            consensus_chars[i] = bases[max_idx]
        end
    end
    
    return DegenerateOlig(consensus_chars)
end

function consensus_degen(msa::AbstractMSA, slack::Real=0.0)
    @assert zero(slack) <= slack <= one(slack) "slack must be between 0.0 and 1.0"
    
    n_positions = length(msa)
    if n_positions == 0
        return DegenerateOlig()
    end
    
    _SORTED_NUC2IUPAC = Dict(
        "A" => 'A',
        "C" => 'C',
        "G" => 'G',
        "T" => 'T',
        "AC" => 'M',
        "AG" => 'R',
        "AT" => 'W',
        "CG" => 'S',
        "CT" => 'Y',
        "GT" => 'K',
        "ACG" => 'V',
        "ACT" => 'H',
        "AGT" => 'D',
        "CGT" => 'B',
        "ACGT" => 'N',
    )
    
    consensus_chars = Vector{Char}(undef, n_positions)
    bases = ('A', 'C', 'G', 'T')
    
    for (i, col) in enumerate(msa)
        max_freq = maximum(col)
        
        if max_freq < slack
            consensus_chars[i] = 'N'
        else
            if slack > zero(slack)
                qualifying_nucs = [bases[j] for j in 1:4 if col[j] >= slack]
                key = String(sort!(qualifying_nucs))
            else
                present_nucs = [bases[j] for j in 1:4 if col[j] > 0.0]
                key = String(sort!(present_nucs))
            end
            
            if isempty(key)
                consensus_chars[i] = 'N'
            else
                consensus_chars[i] = get(_SORTED_NUC2IUPAC, key, 'N')
            end
        end
    end
    
    return DegenerateOlig(consensus_chars)
end



















struct Primer <: AbstractOlig
    seq::String
    gc::Float64
    Tm::Float64
    dg::Float64
    interval::UnitRange{Int64}
    is_forward::Bool
    n_unique_oligs::Int
    msa::AbstractMSA
end

const EMPTY_PRIMER = Primer("", NaN, NaN, NaN, 0:0, false, 1, MSA())
Primer() = EMPTY_PRIMER

isorphan(p::Primer) = isempty(p.msa)

function Primer(
    msa::AbstractMSA,
    is_forward::Bool,
    tail::Integer=NONDEGEN_PRIMER_TAIL,
    head_slack::Real=HEAD_SLACK,
    tail_slack::Real=TAIL_SLACK,
    gc_range::UnitRange{Int}=GC_RANGE,
    tm_range::UnitRange{Int}=TM_RANGE,
    delta_g::Real=MINIMAL_DELTA_G,
    max_olig_vars::Integer=MAX_NONDEG_OLIG_VARIANTS
    )
    @assert zero(head_slack) <= head_slack one(head_slack) "Head slack must be within [0,1]"
    @assert zero(tail_slack) <= tail_slack one(tail_slack) "Tail slack must be within [0,1]"
    interval = 1:length(msa)
    @assert 0<=tail<=length(interval) "Primer tail length must be within [0, $(length(interval))]"
    head_subrange, tail_subrange = is_forward ?
        (range(interval.start, interval.stop-tail), range(interval.stop-tail+1, interval.stop)) : 
        (range(interval.start+tail, interval.stop), range(interval.start, interval.start+tail-1))
    head_seq = consensus_degen(msa[head_subrange], head_slack)
    tail_seq = consensus_major(msa[tail_subrange], tail_slack)
    primer_olig = is_forward ? head_seq*tail_seq : revcomp(tail_seq*head_seq)
    (isempty(primer_olig) || n_unique_oligs(primer_olig) > max_olig_vars) && return Primer()
    
    gc = gc_content(primer_olig)
    (100gc < gc_range.start || 100gc > gc_range.stop)  && return Primer()

    expanded = list_nondegen(primer_olig)
    Tm = sum(tm, expanded) / length(expanded)
    (Tm < tm_range.start || Tm > tm_range.stop)  && return Primer()

    ΔG = minimum(dg, expanded)
    ΔG < delta_g  && return Primer()

    return Primer(
        primer_olig.seq,
        gc,
        Tm,
        ΔG,
        getrange(msa),
        is_forward,
        n_unique_oligs(primer_olig),
        msa
    )
end

function Primer(
    msa::AbstractMSA,
    interval::UnitRange{Int},
    is_forward::Bool,
    tail::Integer=NONDEGEN_PRIMER_TAIL,
    head_slack::Real=HEAD_SLACK,
    tail_slack::Real=TAIL_SLACK,
    gc_range::UnitRange{Int}=GC_RANGE,
    tm_range::UnitRange{Int}=TM_RANGE,
    delta_g::Real=MINIMAL_DELTA_G,
    max_olig_vars::Integer=MAX_NONDEG_OLIG_VARIANTS
    )
    return Primer(
        msa[interval],
        is_forward,
        tail,
        head_slack,
        tail_slack,
        gc_range,
        tm_range,
        delta_g,
        max_olig_vars
    )
end


function Primer(
    primer_seq::AbstractString,
    is_forward::Bool,
    gc_range::UnitRange{Int}=GC_RANGE,
    tm_range::UnitRange{Int}=TM_RANGE,
    delta_g::Real=MINIMAL_DELTA_G,
    max_olig_vars::Integer=MAX_NONDEG_OLIG_VARIANTS
    )
    primer_olig = DegenerateOlig(primer_seq)
    gc = gc_content(primer_olig)
    (isempty(primer_olig) || n_unique_oligs(primer_olig) > max_olig_vars ) && return new()
    
    gc = gc_content(primer_olig)
    (100gc < gc_range.start || 100gc > gc_range.stop)  && return Primer()

    expanded = list_nondegen(primer_olig)
    Tm = sum(tm, expanded) / length(expanded)
    (Tm < tm_range.start || Tm > tm_range.stop)  && return Primer()

    ΔG = minimum(dg, expanded)
    ΔG < delta_g  && return Primer()
    Primer(
        primer_olig.seq,
        gc,
        Tm,
        ΔG,
        1:length(primer_olig),
        is_forward,
        n_unique_oligs(primer_olig),
        MSA()
    )
end


function Base.show(io::IO, primer::Primer)
    if isempty(primer)
        print(io, "Primer()")
        return
    end
    if primer.is_forward
        println(io, "$(primer.seq): forward primer.")
    else
        println(io, "$(primer.seq): reverse primer.")
    end
    print(io, "MSA interval=$(primer.interval), ")
    print(io, "GC=$(round(primer.gc, digits=2)), ")
    print(io, "ΔG=$(round(primer.dg, digits=2)), ")
    println(io, "Tm=$(round(primer.Tm, digits=1))")
end



function _allranges_lazy(msa_length::Integer, primer_length::Integer)
    candidates = (
        i:(i + primer_length - 1)
        for i in 1:(msa_length - primer_length + 1)
        if primer_length <= msa_length
    )
    return candidates
end

function _allranges_lazy(msa_length::Integer, primer_length_range::UnitRange{<:Integer})
    candidates = Iterators.flatten(
        _allranges_lazy(msa_length, primer_length)
        for primer_length in primer_length_range
    )
    return candidates
end

function _allranges(msa_length, primer_length_or_range)
    candidates_generator = _allranges_lazy(msa_length, primer_length_or_range)
    candidates = collect(candidates_generator)
    return candidates
end

function _allranges(msa_length::Integer, primer_length_range::UnitRange{<:Integer})
    all_candidates = UnitRange{Int}[]

    for primer_length in primer_length_range
        candidates = _allranges(msa_length, primer_length)
        append!(all_candidates, candidates)
    end

    return all_candidates
end

function _msa_valid_ranges(
    MSA::AbstractMSA,
    len_range::UnitRange{<:Integer},
    tail::Integer,
    tail_degen::Integer,
    head_degen::Integer,
    maximal_gap_proportion::Real,
    head_degen_slack::Real,
    tail_degen_slack::Real
    )
    degen_mtr = get_base_count(MSA)
    msa_len = length(MSA)
    
    fwd_candidates = UnitRange{Int}[]
    rev_candidates = UnitRange{Int}[]
    
    function __push_candidate!(_candidates, _degen_mtr, _range, _head, _tail_degen, _head_degen, _head_degen_slack, _tail_degen_slack)

        _head_range = _range.start:(_range.start+_head-1)
        _tail_range = (_head_range.stop+1):_range.stop
        
        _head_deg = count(>(1), count.(>(_head_degen_slack), eachcol(_degen_mtr[:, _head_range])))
        _tail_deg = count(>(1), count.(>(_tail_degen_slack), eachcol(_degen_mtr[:, _tail_range])))
        
        if _head_deg <= _head_degen && _tail_deg <= _tail_degen
            push!(_candidates, _range)
        end
    end

    for range in _allranges(msa_len, len_range)
        if any(sum(view(degen_mtr, :, range); dims=1) .< maximal_gap_proportion)
            continue
        end
        tail = min(tail, length(range))
        head = length(range) - tail

        __push_candidate!(fwd_candidates, degen_mtr, range, head, tail_degen, head_degen, head_degen_slack, tail_degen_slack)
        __push_candidate!(rev_candidates, degen_mtr, range, head, tail_degen, head_degen, head_degen_slack, tail_degen_slack)
    end
    
    return fwd_candidates, rev_candidates
end

function construct_primers(
    MSA::AbstractMSA;
    len_range=PRIMER_LEN_RANGE,
    tail=NONDEGEN_PRIMER_TAIL,
    tail_degen=TAIL_DEGEN,
    head_degen=HEAD_DEGEN,
    gc_range=GC_RANGE,
    tm_range=TM_RANGE,
    delta_g=MINIMAL_DELTA_G,
    maximal_gap_proportion=MAXIMAL_GAP_PROPORTION,
    head_slack=HEAD_SLACK,
    tail_slack=TAIL_SLACK,
    max_olig_variants=MAX_NONDEG_OLIG_VARIANTS
    )

    fwd_ranges, rev_ranges = _msa_valid_ranges(
        MSA,
        len_range,
        tail,
        tail_degen,
        head_degen,
        maximal_gap_proportion,
        head_slack,
        tail_slack
    )

    # gap_mask = msa_gap_mask(MSA)
    # fwd_oligos = _all_oligos_from_strand(MSA, fwd_ranges, gap_mask, false, tail, head_degen_slack, tail_degen_slack)
    # rev_oligos = _all_oligos_from_strand(MSA, rev_ranges, gap_mask, true, tail, head_degen_slack, tail_degen_slack)
    fwd_primers = [
        Primer(MSA,interval,true,tail,head_slack,tail_slack,gc_range,tm_range,delta_g,max_olig_variants)
        for interval in fwd_ranges
    ] |> filter(!isempty)

    rev_primers = [
        Primer(MSA,interval,false,tail,head_slack,tail_slack,gc_range,tm_range,delta_g,max_olig_variants)
        for interval in rev_ranges
    ] |> filter(!isempty)

    return fwd_primers, rev_primers
end


function olig_stats(degen_olig::AbstractOlig)
    gc = gc_content(degen_olig)
    expanded = list_nondegen(degen_olig)
    Tm = sum(tm, expanded) / length(expanded)
    ΔG = minimum(dg, expanded)

    return (Tm=Tm, gc=gc, dg=ΔG)
end

function olig_stats(olig::Olig)
    gc = gc_content(olig)
    Tm = tm(olig)
    ΔG = dg(olig)

    return (Tm=Tm, gc=gc, dg=ΔG)
end


