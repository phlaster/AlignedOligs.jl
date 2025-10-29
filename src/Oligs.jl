abstract type AbstractOlig <: AbstractString end

struct Olig <: AbstractOlig
    seq::String
    description::String

    function Olig(seq::AbstractString, description::Union{AbstractString, Integer}="")
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

Olig() = EMPTY_OLIG

Olig(chars::Vector{Char}, descr::AbstractString = "") = Olig(String(chars), descr)

function Olig(olig::AbstractOlig)
    if n_deg_pos(olig) > 0
        throw(InexactError(:Olig, Olig, olig))
    end
    return Olig(String(olig), description(olig))
end

struct DegenerateOlig <: AbstractOlig
    seq::String
    n_deg_pos::Int
    n_unique_oligs::BigInt
    description::String
end

DegenerateOlig() = EMPTY_OLIG

function DegenerateOlig(seq::AbstractString, descr::Union{AbstractString, Integer} = "")
    isempty(seq) && return DegenerateOlig()

    seq = uppercase(seq)
    seq_chars = Set(seq)
    if !issubset(seq_chars, ALL_BASES)
        error("DegenerateOlig contains unallowed characters: $(join(setdiff(seq_chars, ALL_BASES), ", "))")
    end
    
    n_degenerate = count(char -> char in DEGEN_BASES, seq)
    n_possible = reduce(*, (IUPAC_COUNTS[char] for char in seq), init=BigInt(1))
    return DegenerateOlig(seq, n_degenerate, n_possible, string(descr))
end

DegenerateOlig(chars::Vector{Char}, descr::AbstractString = "") = DegenerateOlig(String(chars), descr)
DegenerateOlig(olig::AbstractOlig) = DegenerateOlig(String(olig), n_deg_pos(olig), n_unique_oligs(olig), description(olig))

struct GappedOlig{T<:Union{Olig,DegenerateOlig}} <: AbstractOlig
    parent::T
    gaps::Vector{Pair{Int}}  # sorted, non-overlapping: start position (1-based in ungapped) => length
    total_length::Int

    function GappedOlig(parent::T, gaps::Vector{Pair{Int}}) where T <: Union{Olig,DegenerateOlig}
        parent_len = length(parent)
        if isempty(gaps)
            return new{T}(parent, gaps, parent_len)
        end

        sorted_gaps = sort(gaps; by=first)
        
        for i in 1:length(sorted_gaps)
            start, len = sorted_gaps[i]
            if start < 1 || start > parent_len + 1
                throw(ArgumentError("Gap start position $start must be between 1 and $(parent_len + 1)"))
            end
            if len <= 0
                throw(ArgumentError("Gap length must be positive, got $len"))
            end
            if i > 1 && sorted_gaps[i-1].first == start
                throw(ArgumentError("Duplicate gap start position: $start"))
            end
        end

        total_length = parent_len + sum(len for (_, len) in sorted_gaps)
        return new{T}(parent, sorted_gaps, total_length)
    end
end

GappedOlig() = EMPTY_OLIG

function GappedOlig(seq::AbstractString, descr::AbstractString = "")
    parent_seq = filter(c -> c != '-', seq)
    underlying_olig = try
        Olig(parent_seq, descr)
    catch err
        if err isa ErrorException && startswith(err.msg, "Olig contains unallowed characters")
            DegenerateOlig(parent_seq, descr)
        else
            rethrow(err)
        end
    end

    gaps = Pair{Int}[]
    parent_pos = 0
    i = 1
    n = length(seq)
    
    while i <= n
        if seq[i] != '-'
            parent_pos += 1
            i += 1
        else
            gap_start = parent_pos + 1
            len = 0
            while i <= n && seq[i] == '-'
                len += 1
                i += 1
            end
            push!(gaps, gap_start => len)
        end
    end
    return GappedOlig(underlying_olig, gaps)
end

struct OligView{T<:AbstractOlig} <: AbstractOlig
    parent::T
    range::UnitRange{Int}
end

struct NonDegenIterator{T<:AbstractOlig}
    olig::T
    n_variants::Integer
end

function Base.getindex(olig::AbstractOlig, r::UnitRange{Int})
    @boundscheck checkbounds(1:lastindex(olig), r)
    if olig isa OligView
        adjusted_start = olig_range(olig).start + r.start - 1
        adjusted_range = adjusted_start:(adjusted_start + length(r) - 1)
    else
        adjusted_range = r
    end
    P = parent(olig)
    return OligView(P, adjusted_range)
end
function Base.getindex(olig::T, i::Int) where {T<:AbstractOlig}
    @boundscheck checkbounds(1:lastindex(olig), i)
    String(olig)[i]
end

function Base.getindex(go::GappedOlig, i::Int)
    @boundscheck checkbounds(1:length(go), i)
    cum_gaps = 0
    for (start, len) in go.gaps
        gapped_start = start + cum_gaps
        gapped_end = gapped_start + len - 1
        if gapped_end < i
            cum_gaps += len
            continue
        end
        if gapped_start <= i <= gapped_end
            return '-'
        end
    end
    ungapped_i = i - cum_gaps
    @boundscheck checkbounds(1:length(go.parent), ungapped_i)
    return go.parent[ungapped_i]
end
function Base.getindex(go::GappedOlig, r::UnitRange{Int})
    @boundscheck checkbounds(1:length(go), r)
    return OligView(go, r)
end
function Base.getindex(ov::OligView, i::Int)
    @boundscheck checkbounds(1:length(ov), i)
    return ov.parent[ov.range.start + i - 1]
end

Base.String(olig::Olig) = olig.seq
Base.String(olig::DegenerateOlig) = olig.seq
Base.String(go::GappedOlig) = _build_gapped_string(go)
Base.String(ov::OligView) = String(parent(ov))[olig_range(ov)]

Base.:(==)(o::AbstractOlig, s::SubString{<:Base.AnnotatedString}) = String(o) == String(s)
Base.:(==)(s::SubString{<:Base.AnnotatedString}, o::AbstractOlig) = o == s
Base.:(==)(o::AbstractOlig, s::Base.AnnotatedString) = String(o) == String(s)
Base.:(==)(s::Base.AnnotatedString, o::AbstractOlig) = o == s
Base.:(==)(a::GappedOlig, b::GappedOlig) = (a.parent == b.parent) && (a.gaps == b.gaps) && (a.total_length == b.total_length)

function _build_gapped_string(go::GappedOlig)
    parent_len = length(parent(go))
    buffer = Vector{Char}(undef, length(go))
    ungapped_pos = 1 
    gap_idx = 1
    buffer_idx = 1
    
    while buffer_idx <= go.total_length && ungapped_pos <= parent_len
        if gap_idx <= length(go.gaps) && ungapped_pos == go.gaps[gap_idx].first
            start, len = go.gaps[gap_idx]
            for i in buffer_idx:(buffer_idx + len - 1)
                buffer[i] = '-'
            end
            buffer_idx += len
            gap_idx += 1
            continue
        end
        buffer[buffer_idx] = parent(go).seq[ungapped_pos]
        ungapped_pos += 1
        buffer_idx += 1
    end
    
    while buffer_idx <= go.total_length && gap_idx <= length(go.gaps)
        start, len = go.gaps[gap_idx]
        for i in buffer_idx:(buffer_idx + len - 1)
            buffer[i] = '-'
        end
        buffer_idx += len
        gap_idx += 1
    end
    
    if ungapped_pos <= parent_len || buffer_idx <= go.total_length
        throw(ArgumentError("Mismatch in sequence construction: ungapped_pos=$ungapped_pos (expected $(parent_len + 1)), buffer_idx=$buffer_idx (expected $(go.total_length + 1))"))
    end
    return String(buffer)
end

Base.parent(olig::AbstractOlig) = olig
Base.parent(olig::OligView) = olig.parent
Base.parent(olig::GappedOlig) = olig.parent
Base.parent(ndi::NonDegenIterator) = ndi.olig

Base.ncodeunits(olig::AbstractOlig) = length(olig)
Base.codeunit(olig::Olig, i::Integer) = codeunit(String(olig), i)
Base.codeunit(olig::DegenerateOlig, i::Integer) = codeunit(String(olig), i)
Base.codeunit(go::GappedOlig, i::Integer) = UInt8(_char_at(go, i))
Base.codeunit(ov::OligView, i::Integer) = codeunit(parent(ov), olig_range(ov).start + i - 1)

function _char_at(go::GappedOlig, k::Integer)
    @boundscheck 1 <= k <= go.total_length || throw(BoundsError(go, k))
    parent_len = length(parent(go))
    cum_gaps = 0
    ungapped_pos = 0 
    gapped_pos = 0
    
    for (start, len) in go.gaps
        positions_before_gap = start - ungapped_pos
        gapped_end = gapped_pos + positions_before_gap
        if gapped_pos < k <= gapped_end
            return String(parent(go))[ungapped_pos + (k - gapped_pos)]
        end
        gapped_pos = gapped_end
        if gapped_pos < k <= gapped_pos + len
            return '-'
        end
        gapped_pos += len
        cum_gaps += len
        ungapped_pos = start
    end
    
    if gapped_pos < k <= gapped_pos + (parent_len - ungapped_pos)
        return String(parent(go))[ungapped_pos + (k - gapped_pos)]
    end
    throw(BoundsError(go, k))
end

Base.isvalid(olig::AbstractOlig, i::Int) = 1 <= i <= length(olig)

Base.length(olig::Olig) = length(String(olig))
Base.length(olig::DegenerateOlig) = length(String(olig))
Base.length(go::GappedOlig) = go.total_length
Base.length(ov::OligView) = length(olig_range(ov))
Base.isempty(olig::AbstractOlig) = length(olig) == 0
Base.lastindex(olig::AbstractOlig) = length(olig)

olig_range(olig::AbstractOlig) = 1:length(olig)
olig_range(ov::OligView) = ov.range

description(olig::AbstractOlig) = parent(olig).description
description(ov::OligView{<:GappedOlig}) = parent(parent(ov)).description

hasgaps(::AbstractOlig) = false
hasgaps(ov::OligView{Olig}) = false
hasgaps(ov::OligView{DegenerateOlig}) = false
hasgaps(go::GappedOlig) = !isempty(go.gaps)
hasgaps(ov::OligView) = any(c == '-' for c in ov)


n_unique_oligs(::AbstractOlig) = BigInt(1)
n_unique_oligs(d::DegenerateOlig) = d.n_unique_oligs
n_unique_oligs(ov::OligView) = reduce(*, (IUPAC_COUNTS[base] for base in ov), init=BigInt(1))
n_unique_oligs(go::GappedOlig) = n_unique_oligs(parent(go))

n_deg_pos(::AbstractOlig) = 0
n_deg_pos(d::DegenerateOlig) = d.n_deg_pos
n_deg_pos(ov::OligView) = count(char -> char in DEGEN_BASES, ov)
n_deg_pos(go::GappedOlig) = n_deg_pos(parent(go))

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

function Base.iterate(iter::NonDegenIterator)
    olig = parent(iter)
    length(olig) == 0 && return nothing
    
    options = [IUPAC_B2V[c] for c in String(olig)]
    lens = length.(options)
    indices = ones(Int, length(olig))
    buffer = Vector{Char}(undef, length(olig))
    
    for j in 1:length(olig)
        buffer[j] = options[j][indices[j]]
    end
    
    state = (indices, options, lens, buffer, length(olig))
    descr = string("Non-degen sample from: ", description(olig))
    return Olig(String(buffer), descr), state
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
    descr = string("Non-degen sample from: ", description(parent(iter)))
    return Olig(String(buffer), descr), (indices, options, lens, buffer, n)
end

Base.length(iter::NonDegenIterator) = iter.n_variants
Base.eltype(::Type{<:NonDegenIterator}) = Olig

nondegens(olig::Olig) = (olig,)
nondegens(go::GappedOlig{Olig}) = (go,)
nondegens(olig::DegenerateOlig) = NonDegenIterator(olig, n_unique_oligs(olig))
function nondegens(go::GappedOlig{DegenerateOlig})
    (GappedOlig(olig, go.gaps) for olig in nondegens(go.parent))
end
function nondegens(ov::OligView)
    str = String(ov)
    descr = description(ov)
    if hasgaps(ov)
        go = GappedOlig(str, descr)
        return nondegens(go)
    else
        deg_olig = DegenerateOlig(str, descr)
        return nondegens(deg_olig)
    end
end

Base.rand(rng::AbstractRNG, olig::Olig) = olig
function Base.rand(rng::AbstractRNG, olig::DegenerateOlig)
    isempty(olig) && return EMPTY_OLIG
    buffer = Vector{Char}(undef, length(olig))
    @inbounds for (i, c) in enumerate(olig)
        options = IUPAC_B2V[c]
        buffer[i] = rand(rng, options)
    end
    descr = string("Random non-degen sample from: ", description(olig))
    return Olig(String(buffer), descr)
end
function Base.rand(rng::AbstractRNG, ov::OligView)
    parent_olig = rand(rng, parent(ov))
    descr = string("Random non-degen sample from $(olig_range(ov)) OligView of: ", description(ov))
    return Olig(String(parent_olig)[olig_range(ov)], descr)
end
Base.rand(::AbstractRNG, go::GappedOlig) = go
Base.rand(olig::AbstractOlig) = rand(Random.GLOBAL_RNG, olig)

function Base.iterate(go::GappedOlig)
    length(go) == 0 && return nothing
    return iterate(go, (1, 1, 1, 0))
end
function Base.iterate(go::GappedOlig, state::Tuple{Int, Int, Int, Int})
    pos, ungapped_pos, gap_idx, remaining = state
    if pos > length(go)
        return nothing
    end
    if remaining > 0
        char = '-'
        new_remaining = remaining - 1
        new_ungapped = ungapped_pos
        new_gap_idx = new_remaining == 0 ? gap_idx + 1 : gap_idx
    elseif gap_idx <= length(go.gaps) && ungapped_pos == go.gaps[gap_idx].first
        start, len = go.gaps[gap_idx]
        char = '-'
        new_remaining = len - 1
        new_ungapped = ungapped_pos
        new_gap_idx = new_remaining == 0 ? gap_idx + 1 : gap_idx
    else
        if ungapped_pos > length(go.parent)
            throw(BoundsError(go.parent, ungapped_pos))
        end
        char = go.parent[ungapped_pos]
        new_ungapped = ungapped_pos + 1
        new_remaining = 0
        new_gap_idx = gap_idx
    end
    return char, (pos + 1, new_ungapped, new_gap_idx, new_remaining)
end
Base.iterate(olig::AbstractOlig) = length(olig) == 0 ? nothing : (olig[1], 2)
function Base.iterate(olig::AbstractOlig, state::Int)
    if state > length(olig)
        return nothing
    end
    return olig[state], state + 1
end
