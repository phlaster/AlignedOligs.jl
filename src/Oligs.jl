abstract type AbstractOlig <: AbstractString end

struct Olig <: AbstractOlig
    seq::String
    description::String

    function Olig(seq::AbstractString, description = "")
        isempty(seq) && return new("", "")
        seq = uppercase(seq)
        seq_chars = Set(seq)
        if !issubset(seq_chars, NON_DEGEN_BASES)
            error("Olig contains unallowed characters: $(join(setdiff(seq_chars, NON_DEGEN_BASES), ", "))")
        end
        return new(seq, string(description))
    end
end

const EMPTY_OLIG = Olig("", "")

struct DegenerateOlig <: AbstractOlig
    seq::String
    n_deg_pos::Int
    n_unique_oligs::BigInt
    description::String
end

struct GappedOlig <: AbstractOlig
    olig::Olig
    gaps::Vector{Pair{Int, Int}}  # sorted, non-overlapping: start position (1-based in gapped) => length
    total_length::Int
    description::String

    function GappedOlig(olig::Olig, gaps::Vector{Pair{Int, Int}}, descr = "")
        new_descr = isempty(descr) ? description(olig) : descr
        isempty(gaps) && return new(olig, gaps, length(olig), new_descr)
        sort!(gaps, by=first)  # Ensure sorted by start
        prev_end = 0
        total_gaps = 0
        for (start, len) in gaps
            start >= 1 || throw(ArgumentError("Gap start must be at least 1"))
            len > 0 || throw(ArgumentError("Gap length must be positive"))
            start > prev_end || throw(ArgumentError("Gaps must be non-overlapping"))
            prev_end = start + len - 1
            total_gaps += len
        end
        total_len = length(olig) + total_gaps
        prev_end <= total_len || throw(ArgumentError("Gaps exceed sequence length"))
        new(olig, gaps, total_len, new_descr)
    end
end

function GappedOlig(seq::AbstractString, descr::AbstractString="")
    gaps = Pair{Int, Int}[]
    i = 1
    n = length(seq)
    while i <= n
        if seq[i] == '-'
            start = i
            len = 1
            i += 1
            while i <= n && seq[i] == '-'
                len += 1
                i += 1
            end
            push!(gaps, start => len)
        else
            i += 1
        end
    end
    
    non_gapped_seq = replace(seq, "-" => "")
    underlying_olig = Olig(non_gapped_seq)
    
    return GappedOlig(underlying_olig, gaps, descr)
end

struct OligView{T<:AbstractOlig} <: AbstractOlig
    parent::T
    range::UnitRange{Int}
end
    
Base.String(olig::Olig) = olig.seq
Base.String(olig::DegenerateOlig) = olig.seq
Base.String(go::GappedOlig) = _build_gapped_string(go)
Base.String(ov::OligView) = String(parent(ov))[ov.range]

Base.:(==)(o::AbstractOlig, s::SubString{<:Base.AnnotatedString}) = String(o) == String(s)
Base.:(==)(s::SubString{<:Base.AnnotatedString}, o::AbstractOlig) = o == s
Base.:(==)(o::AbstractOlig, s::Base.AnnotatedString) = String(o) == String(s)
Base.:(==)(s::Base.AnnotatedString, o::AbstractOlig) = o == s


function _build_gapped_string(go::GappedOlig)
    buffer = Char[]
    ungapped_pos = 1
    gap_idx = 1
    for pos in 1:go.total_length
        if gap_idx <= length(go.gaps)
            start, len = go.gaps[gap_idx]
            if start <= pos < start + len
                push!(buffer, '-')
                continue
            elseif pos == start + len
                gap_idx += 1
            end
        end
        push!(buffer, go.olig.seq[ungapped_pos])
        ungapped_pos += 1
    end
    return String(buffer)
end

parent(olig::AbstractOlig) = olig
parent(olig::OligView) = olig.parent

Base.ncodeunits(olig::AbstractOlig) = length(olig)
Base.codeunit(olig::Olig, i::Integer) = codeunit(String(olig), i)
Base.codeunit(olig::DegenerateOlig, i::Integer) = codeunit(String(olig), i)
Base.codeunit(go::GappedOlig, i::Integer) = UInt8(_char_at(go, i))
Base.codeunit(ov::OligView, i::Integer) = codeunit(parent(ov), ov.range.start + i - 1)

function _char_at(go::GappedOlig, k::Integer)
    @boundscheck 1 <= k <= go.total_length || throw(BoundsError(go, k))
    cum_gaps = 0
    for (start, len) in go.gaps
        if start + len <= k
            cum_gaps += len
        elseif start <= k
            return '-'
        end
    end
    ungapped_pos = k - cum_gaps
    return go.olig.seq[ungapped_pos]
end

Base.isvalid(olig::AbstractOlig, i::Int) = 1 <= i <= length(olig)

Base.lastindex(olig::AbstractOlig) = length(olig)
Base.length(olig::Olig) = length(String(olig))
Base.length(olig::DegenerateOlig) = length(String(olig))
Base.length(go::GappedOlig) = go.total_length
Base.length(ov::OligView) = length(ov.range)
Base.isempty(olig::AbstractOlig) = length(olig) == 0

description(olig::AbstractOlig) = parent(olig).description

Olig(chars::Vector{Char}, descr::AbstractString = "") = Olig(String(chars), descr)

function Olig(olig::AbstractOlig)
    if n_deg_pos(olig) > 0
        throw(InexactError(:Olig, Olig, olig))
    end
    return Olig(String(olig), description(olig))
end

Olig() = EMPTY_OLIG

function DegenerateOlig(seq::AbstractString, descr = "")
    isempty(seq) && return DegenerateOlig()

    seq = uppercase(seq)
    seq_chars = Set(seq)
    if !issubset(seq_chars, ALL_BASES)
        error("DegenerateOlig contains unallowed characters: $(join(setdiff(seq_chars, ALL_BASES), ", "))")
    end
    
    n_degenerate = count(char -> char in DEGEN_BASES, seq)
    n_possible = reduce(*, IUPAC_COUNTS[char] for char in seq, init=BigInt(1))
    return DegenerateOlig(seq, n_degenerate, n_possible, string(descr))
end

DegenerateOlig() = EMPTY_OLIG

DegenerateOlig(chars::Vector{Char}, descr::AbstractString = "") = DegenerateOlig(String(chars), descr)
DegenerateOlig(olig::AbstractOlig) = DegenerateOlig(String(olig), n_deg_pos(olig), n_unique_oligs(olig), description(olig))

function Base.getindex(olig::AbstractOlig, r::UnitRange{Int})
    @boundscheck checkbounds(1:lastindex(olig), r)
    if olig isa OligView
        adjusted_start = olig.range.start + r.start - 1
        adjusted_range = adjusted_start:(adjusted_start + length(r) - 1)
    else
        adjusted_range = r
    end
    P = parent(olig)
    return OligView{typeof(P)}(P, adjusted_range)
end

function Base.getindex(olig::T, i::Int) where {T<:AbstractOlig}
    @boundscheck checkbounds(1:lastindex(olig), i)
    String(olig)[i]
end

hasgaps(::AbstractOlig) = false
hasgaps(go::GappedOlig) = !isempty(go.gaps)
hasgaps(ov::OligView) = hasgaps(parent(ov))

n_unique_oligs(::AbstractOlig) = BigInt(1)
n_unique_oligs(d::DegenerateOlig) = d.n_unique_oligs
n_unique_oligs(ov::OligView) = hasgaps(ov) ? BigInt(1) : reduce(*, (IUPAC_COUNTS[base] for base in ov), init=BigInt(1))

n_deg_pos(::AbstractOlig) = 0
n_deg_pos(d::DegenerateOlig) = d.n_deg_pos
n_deg_pos(ov::OligView) = hasgaps(ov) ? 0 : count(char -> char in DEGEN_BASES, ov)

function Base.:*(o1::AbstractOlig, o2::AbstractOlig)
    if hasgaps(o1) || hasgaps(o2)
        throw(ErrorException("Concatenation not implemented for gapped oligs"))
    end
    if n_deg_pos(o1) > 0 || n_deg_pos(o2) > 0
        return DegenerateOlig(String(o1) * String(o2), n_deg_pos(o1)+n_deg_pos(o2), n_unique_oligs(o1)*n_unique_oligs(o2), "concat")
    else
        return Olig(String(o1) * String(o2), "concat")
    end
end

Base.convert(::Olig, o::AbstractOlig) = Olig(o)
Base.convert(::Type{DegenerateOlig}, o::AbstractOlig) = DegenerateOlig(o)
Base.promote_rule(::Type{Olig}, ::Type{DegenerateOlig}) = DegenerateOlig
Base.promote_rule(::Type{Olig}, ::Type{GappedOlig}) = GappedOlig
Base.promote_rule(::Type{DegenerateOlig}, ::Type{GappedOlig}) = DegenerateOlig

# function Base.show(io::IO, olig::AbstractOlig)
#     color = get(io, :color, false)
#     yellow = color ? "\e[33m" : ""
#     reset = color ? "\e[0m" : ""

#     seq = String(olig)
#     max_width = 20
#     if length(seq) > max_width
#         seq_display = seq[1:min(max_width-3, end)] * "..."
#         seq_info = " $yellow$(length(seq)) nt$reset"
#     else
#         seq_display = seq
#         seq_info = ""
#     end
    
#     if olig isa DegenerateOlig
#         deg_str = n_deg_pos(olig) == 0 ? "no degenerate positions" : "$(n_deg_pos(olig)) degenerate positions"
#         variants = n_unique_oligs(olig)
#         variants_str = variants > 10_000 ? ">10k" : string(variants)
        
#         println(io, "DegenerateOlig: $seq_display$seq_info")
#         println(io, "  • $deg_str, $variants_str unique variants")
#     elseif olig isa GappedOlig
#         gap_count = sum(last.(olig.gaps))
#         gap_str = gap_count == 0 ? "no gaps" : "$gap_count gap positions"
#         println(io, "GappedOlig: $seq_display$seq_info")
#         println(io, "  • $gap_str")
#     elseif olig isa OligView
#         println(io, "OligView: $seq_display$seq_info")
#     else
#         println(io, "Olig: $seq_display$seq_info")
#     end
    
#     descr = description(olig)
#     if !isempty(descr)
#         desc_info = length(descr) > 60 ? " $yellow$(length(descr)) bytes total$reset" : ""
#         desc_display = length(descr) > 42 ? descr[1:min(42, end)] * "..." : descr
#         print(io, "  • $desc_display$desc_info")
#     end
# end

struct NonDegenIterator{T<:AbstractOlig}
    olig::T
    n_variants::Integer
end

parent(ndi::NonDegenIterator) = ndi.olig

function Base.iterate(iter::NonDegenIterator)
    olig = iter.olig
    length(olig) == 0 && return nothing
    
    options = [IUPAC_B2V[c] for c in String(olig)]
    lens = length.(options)
    indices = ones(Int, length(olig))
    buffer = Vector{Char}(undef, length(olig))
    
    for j in 1:length(olig)
        buffer[j] = options[j][indices[j]]
    end
    
    state = (indices, options, lens, buffer, length(olig))
    return Olig(String(buffer), description(olig)), state
end

function Base.iterate(iter::NonDegenIterator, state)
    indices, options, lens, buffer, n = state
    n == 0 && return nothing

    pos = n
    @inbounds while pos > 0
        indices[pos] += 1
        indices[pos] <= lens[pos] && break
        indices[pos] = 1
        pos -= 1
    end
    pos == 0 && return nothing

    @inbounds @simd for j in 1:n
        buffer[j] = options[j][indices[j]]
    end

    return Olig(String(buffer), description(parent(iter))), (indices, options, lens, buffer, n)
end

Base.length(iter::NonDegenIterator) = iter.n_variants
Base.eltype(::Type{<:NonDegenIterator}) = Olig

nondegens(olig::Olig) = (olig,)
nondegens(go::GappedOlig) = (go,)
nondegens(olig::DegenerateOlig) = NonDegenIterator(olig, n_unique_oligs(olig))
nondegens(ov::OligView) = NonDegenIterator(DegenerateOlig(String(ov), description(ov)), n_unique_oligs(ov))

Base.rand(rng::AbstractRNG, olig::Olig) = olig

function Base.rand(rng::AbstractRNG, olig::DegenerateOlig)
    isempty(olig) && return EMPTY_OLIG
    buffer = Vector{Char}(undef, length(olig))
    @inbounds for (i, c) in enumerate(olig)
        options = IUPAC_B2V[c]
        buffer[i] = rand(rng, options)
    end
    return Olig(String(buffer), description(olig))
end

function Base.rand(rng::AbstractRNG, ov::OligView)
    parent_olig = rand(rng, parent(ov))
    return Olig(String(parent_olig)[ov.range], description(ov))
end

Base.rand(rng::AbstractRNG, go::GappedOlig) = go

Base.rand(olig::AbstractOlig) = rand(Random.GLOBAL_RNG, olig)

function Base.iterate(go::GappedOlig)
    length(go) == 0 && return nothing
    return iterate(go, (1, 1, 1))
end

function Base.iterate(go::GappedOlig, state::Tuple{Int, Int, Int})
    pos, ungapped_pos, gap_idx = state
    if pos > length(go)
        return nothing
    end
    char = nothing
    new_ungapped = ungapped_pos
    new_gap_idx = gap_idx
    if gap_idx <= length(go.gaps)
        start, len = go.gaps[gap_idx]
        if start <= pos < start + len
            char = '-'
        elseif pos == start + len
            char = go.olig.seq[ungapped_pos]
            new_gap_idx += 1
            new_ungapped += 1
        else
            char = go.olig.seq[ungapped_pos]
            new_ungapped += 1
        end
    else
        char = go.olig.seq[ungapped_pos]
        new_ungapped += 1
    end
    return char, (pos + 1, new_ungapped, new_gap_idx)
end

Base.iterate(olig::AbstractOlig) = length(olig) == 0 ? nothing : (olig[1], 2)

function Base.iterate(olig::AbstractOlig, state::Int)
    if state > length(olig)
        return nothing
    end
    return olig[state], state + 1
end

# SeqFold specifications
SeqFold.revcomp(olig::Olig) = Olig(SeqFold.revcomp(String(olig); table=DNA_COMP_TABLE_DEG))
SeqFold.revcomp(olig::DegenerateOlig) = DegenerateOlig(SeqFold.revcomp(String(olig); table=DNA_COMP_TABLE_DEG))
function SeqFold.revcomp(go::GappedOlig)
    rev_olig = SeqFold.revcomp(go.olig)
    total_len = length(go)
    rev_gaps = Pair{Int, Int}[]
    for (start, len) in go.gaps
        new_start = total_len - start - len + 2
        push!(rev_gaps, new_start => len)
    end
    sort!(rev_gaps, by=first)
    return GappedOlig(rev_olig, rev_gaps)
end
SeqFold.revcomp(ov::OligView) = SeqFold.revcomp(parent(ov))[length(parent(ov)) - ov.range.stop + 1 : length(parent(ov)) - ov.range.start + 1]  # Reverse view

SeqFold.complement(olig::Olig) = Olig(SeqFold.complement(String(olig); table=DNA_COMP_TABLE_DEG))
SeqFold.complement(olig::DegenerateOlig) = DegenerateOlig(SeqFold.complement(String(olig); table=DNA_COMP_TABLE_DEG))
SeqFold.complement(go::GappedOlig) = GappedOlig(SeqFold.complement(go.olig), copy(go.gaps))
SeqFold.complement(ov::OligView) = SeqFold.complement(parent(ov))[ov.range]

function SeqFold.gc_content(olig::AbstractOlig)::Float64
    if hasgaps(olig)
        return SeqFold.gc_content(olig.olig)  # Ignore gaps
    end
    if isempty(olig)
        return NaN
    end
    total_gc = 0.0
    for c in olig
        total_gc += IUPAC_GC_CONTENT[c]
    end
    
    return total_gc / length(olig)
end

function SeqFold.fold(olig::AbstractOlig; temp::Real = 37.0)::Vector{Structure} 
    hasgaps(olig) && throw(ErrorException("Folding not supported for gapped sequences"))
    return SeqFold.fold(String(Olig(olig)); temp=temp)
end

function SeqFold.dg(olig::AbstractOlig; temp::Real = 37.0)::Float64 
    hasgaps(olig) && throw(ErrorException("Free energy calculation not supported for gapped sequences"))
    return SeqFold.dg(String(Olig(olig)); temp=temp)
end

function SeqFold.dg_cache(olig::AbstractOlig; temp::Real = 37.0)::Matrix{Float64}
    if hasgaps(olig)
        throw(ErrorException("Free energy cache not supported for gapped sequences"))
    else
        return SeqFold.dg_cache(String(Olig(olig)); temp=temp)
    end
end

SeqFold.dot_bracket(olig::AbstractOlig, structs::Vector{SeqFold.Structure}) = SeqFold.dot_bracket(String(Olig(olig)), structs)

function SeqFold.tm(olig1::AbstractOlig, olig2::AbstractOlig; conditions=:pcr, conf_int::Real=0.8, max_variants::Int=10000, kwargs...)
    (hasgaps(olig1) || hasgaps(olig2)) && throw(ErrorException("Melting temperature calculation not supported for gapped sequences"))
    if !(0 < conf_int <= 1)
        throw(ArgumentError("conf_int must be in range (0, 1]"))
    end

    i = n_unique_oligs(olig1)
    j = n_unique_oligs(olig2)

    if i == 1 && j == 1
        mean_Tm = SeqFold.tm(String(olig1), String(olig2); conditions=conditions, kwargs...)
        return (
            mean = mean_Tm,
            conf = (mean_Tm, mean_Tm)
        )
    end

    T = zeros(Float64, min(i * j, max_variants))
    if i * j > max_variants
        for k in 1:max_variants
            o1 = rand(olig1)
            o2 = rand(olig2)
            @inbounds T[k] = SeqFold.tm(String(o1), String(o2); conditions=conditions, kwargs...)
        end
    else
        counter = 1
        for o1 in nondegens(olig1)
            for o2 in nondegens(olig2)
                @inbounds T[counter] = SeqFold.tm(String(o1), String(o2); conditions=conditions, kwargs...)
                counter += 1
            end
        end
    end

    mean_Tm = mean(T)
    alpha = (1 - conf_int) / 2
    low = quantile(T, alpha)
    high = quantile(T, 1 - alpha)

    return (
        mean = round(mean_Tm, digits=1),
        conf = (round(low, digits=1), round(high, digits=1))
    )
end

function SeqFold.tm(olig::AbstractOlig; conditions=:pcr, conf_int::Real=0.9, max_variants::Int=10000, kwargs...)
    SeqFold.tm(olig, SeqFold.complement(olig); conditions=conditions, conf_int=conf_int, max_variants=max_variants, kwargs...)
end

function SeqFold.tm_cache(olig1::AbstractOlig, olig2::AbstractOlig; conditions=:pcr, kwargs...)::Matrix{Float64}
    (hasgaps(olig1) || hasgaps(olig2)) && throw(ErrorException("Melting temperature cache not supported for gapped sequences"))
    SeqFold.tm_cache(String(Olig(olig1)), String(Olig(olig2)); conditions=conditions, kwargs...)
end

function SeqFold.tm_cache(olig::AbstractOlig; conditions=:pcr, kwargs...)::Matrix{Float64}
    SeqFold.tm_cache(olig, SeqFold.complement(olig); conditions=conditions, kwargs...)
end

SeqFold.gc_cache(olig::AbstractOlig)::Matrix{Float64} = SeqFold.gc_cache(String(Olig(olig)))