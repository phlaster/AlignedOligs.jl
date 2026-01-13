module Oligs
include("utils.jl")

using Statistics

export AbstractOlig, AbstractDegen, AbstractGapped
export Olig, DegenOlig, GappedOlig, OligView
export olig_range, description
export NonDegenIterator, nondegens
export hasgaps, getgaps, n_deg_pos, n_unique_oligs
export sampleChar, sampleView, sampleNondeg, sample_max_gc, sample_min_gc
export unfolded_proportion

abstract type AbstractOlig <: AbstractString end
abstract type AbstractDegen <: AbstractOlig end
abstract type AbstractGapped <: AbstractDegen end

description(::AbstractString) = ""
struct Olig <: AbstractOlig
    seq::String
    description::String
    function Olig(seq::T, descr::Union{AbstractString,Integer}="") where T <: AbstractString
        descr = (_d = string(descr); isempty(_d)) ? description(seq) : _d
        seq = uppercase(seq)
        seq_chars = Set(seq)
        if !isempty(seq) && !issubset(seq_chars, NON_DEGEN_BASES)
            error("Olig contains unallowed characters: $(join(setdiff(seq_chars, NON_DEGEN_BASES), ", "))")
        end
        new(seq, descr)
    end
end
struct DegenOlig <: AbstractDegen
    seq::String
    n_deg_pos::Int
    n_unique_oligs::BigInt
    description::String

    function DegenOlig(seq::AbstractString, n_deg::Int, n_unique::BigInt, descr::Union{AbstractString,Integer})
        descr = string(descr)
        seq_chars = Set(seq)
        if !isempty(seq) && !issubset(seq_chars, ALL_BASES)
            error("DegenOlig contains unallowed characters: $(join(setdiff(seq_chars, ALL_BASES), ", "))")
        end
        
        actual_deg = count(char -> char in DEGEN_BASES, seq)
        if actual_deg != n_deg
            error("Inconsistent degenerate position count: calculated $actual_deg, given $n_deg")
        end
        
        actual_unique = isempty(seq) ? BigInt(1) : reduce(*, (IUPAC_COUNTS[char] for char in seq), init=BigInt(1))
        if actual_unique != n_unique
            error("Inconsistent unique olig count: calculated $actual_unique, given $n_unique")
        end
        
        new(uppercase(seq), n_deg, n_unique, descr)
    end
end
struct GappedOlig <: AbstractGapped
    parent::DegenOlig
    gaps::Vector{Pair{Int, Int}}
    total_length::Int
    
    function GappedOlig(parent::DegenOlig, gaps::Vector{Pair{Int, Int}}, total_len::Int)
        issorted(gaps, by=first) || error("Gaps ids are not sorted")
        
        calculated_length = length(parent.seq) + sum(x->x.second, gaps, init=0)
        if calculated_length != total_len
            error("Inconsistent total length: calculated $calculated_length, given $total_len")
        end
        
        new(parent, gaps, total_len)
    end
end

const EMPTY_OLIG = Olig("", "")
Olig() = EMPTY_OLIG

const EMPTY_DEGENERATE = DegenOlig("", 0, BigInt(1), "")
DegenOlig() = EMPTY_DEGENERATE

const EMPTY_GAPPED = GappedOlig(DegenOlig(), Pair{Int, Int}[], 0)
GappedOlig() = EMPTY_GAPPED

function DegenOlig(seq::AbstractString, descr::Union{AbstractString,Integer}="")
    descr = string(descr)
    if isempty(seq)
        if isempty(descr)
            return DegenOlig()
        end
        return DegenOlig("", 0, BigInt(1), descr)
    end

    seq = uppercase(seq)
    seq_chars = Set(seq)
    if !issubset(seq_chars, ALL_BASES)
        error("DegenOlig contains unallowed characters: $(join(setdiff(seq_chars, ALL_BASES), ", "))")
    end
    
    n_degenerate = count(char -> char in DEGEN_BASES, seq)
    n_possible = reduce(*, (IUPAC_COUNTS[char] for char in seq), init=BigInt(1))
    return DegenOlig(seq, n_degenerate, n_possible, descr)
end

function GappedOlig(seq::AbstractString, descr::Union{AbstractString,Integer}="")
    descr = string(descr)
    if isempty(seq)
        if isempty(descr)
            return GappedOlig()
        end
        return GappedOlig(DegenOlig("", descr), Pair{Int, Int}[], 0)
    end
    
    parent_seq = filter(!=('-'), seq)
    underlying_olig = DegenOlig(parent_seq, descr)

    gaps = Pair{Int, Int}[]
    parent_pos = 0
    i = 1
    n = length(seq)
    
    @inbounds while i <= n
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
    return GappedOlig(underlying_olig, gaps, n)
end

struct OligView{T<:AbstractOlig} <: AbstractOlig
    parent::T
    range::UnitRange{Int}

    function OligView(olig::T, interval::UnitRange{Int}) where T<:Union{Olig, DegenOlig, GappedOlig}
        @boundscheck checkbounds(olig, interval)
        new{T}(olig, interval)
    end
end

struct NonDegenIterator{T<:AbstractOlig}
    olig::T
    n_variants::Integer
end

Olig(chars::Vector{Char}, descr::Union{AbstractString,Integer}="") = Olig(String(chars), descr)
Olig(olig::Olig) = olig

DegenOlig(chars::Vector{Char}, descr::Union{AbstractString,Integer}="") = DegenOlig(String(chars), descr)
function DegenOlig(olig::AbstractOlig, descr::Union{AbstractString,Integer}="")
    new_descr = isempty(descr) ? description(olig) : descr
    DegenOlig(String(olig), new_descr)
end
DegenOlig(olig::DegenOlig) = olig

GappedOlig(chars::Vector{Char}, descr::Union{AbstractString,Integer}="") = DegenOlig(String(chars), descr)
function GappedOlig(olig::AbstractOlig, descr::Union{AbstractString,Integer}="")
    new_descr = isempty(descr) ? description(olig) : descr
    GappedOlig(String(olig), new_descr)
end
GappedOlig(olig::GappedOlig) = olig

Base.getindex(olig::AbstractOlig, r::UnitRange{Int}) = OligView(parent(olig), r)
Base.getindex(olig::AbstractOlig, i::Int) = getindex(String(olig), i)
function Base.getindex(ov::OligView, r::UnitRange{Int})
    @boundscheck checkbounds(ov, r)
    ov_range_start::Int = olig_range(ov).start
    actual_start = ov_range_start + r.start - 1
    actual_stop = ov_range_start + r.stop - 1
    actual_range = actual_start:actual_stop
    return OligView(parent(ov), actual_range)
end
function Base.getindex(olig::GappedOlig, r::UnitRange{Int})
    @boundscheck begin
        if !(1 <= first(r) <= last(r) <= length(olig))
            throw(BoundsError(olig, r))
        end
    end
    return OligView(olig, r)
end
function Base.getindex(ov::OligView, i::Int)
    adjusted_i = olig_range(ov).start + i - 1
    @boundscheck checkbounds(ov, i)
    return getindex(parent(ov), adjusted_i)
end

Base.String(olig::Olig) = olig.seq
Base.String(olig::DegenOlig) = olig.seq
Base.String(ov::OligView) = String(parent(ov))[olig_range(ov)]
function Base.String(go::GappedOlig)
    parent_len = length(parent(go))
    buffer = Vector{Char}(undef, length(go))
    ungapped_pos = 1 
    gap_idx = 1
    buffer_idx = 1
    gaps = getgaps(go)
    L = length(go)

    @inbounds while buffer_idx <= L && ungapped_pos <= parent_len
        if gap_idx <= length(gaps) && ungapped_pos == gaps[gap_idx].first
            _, len = gaps[gap_idx]
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
    
    @inbounds while buffer_idx <= L && gap_idx <= length(gaps)
        _, len = gaps[gap_idx]
        for i in buffer_idx:(buffer_idx + len - 1)
            buffer[i] = '-'
        end
        buffer_idx += len
        gap_idx += 1
    end
    
    if ungapped_pos <= parent_len || buffer_idx <= L
        throw(ArgumentError("Mismatch in sequence construction: ungapped_pos=$ungapped_pos (expected $(parent_len + 1)), buffer_idx=$buffer_idx (expected $(L + 1))"))
    end
    return String(buffer)
end

Base.:(==)(o::AbstractOlig, s::SubString{<:Base.AnnotatedString}) = String(o) == String(s)
Base.:(==)(s::SubString{<:Base.AnnotatedString}, o::AbstractOlig) = o == s
Base.:(==)(o::AbstractOlig, s::Base.AnnotatedString) = String(o) == String(s)
Base.:(==)(s::Base.AnnotatedString, o::AbstractOlig) = o == s
Base.:(==)(a::GappedOlig, b::GappedOlig) = (parent(a) == parent(b)) && (getgaps(a) == getgaps(b)) && (length(a) == length(b))

Base.parent(olig::AbstractOlig) = olig
Base.parent(olig::OligView) = olig.parent
Base.parent(olig::GappedOlig) = olig.parent
Base.parent(ndi::NonDegenIterator) = ndi.olig

Base.ncodeunits(olig::AbstractOlig) = length(olig)

Base.codeunit(::AbstractOlig) = UInt8
Base.codeunit(olig::Olig, i::Integer) = codeunit(String(olig), i)
Base.codeunit(olig::DegenOlig, i::Integer) = codeunit(String(olig), i)
Base.codeunit(ov::OligView, i::Integer) = codeunit(parent(ov), olig_range(ov).start + i - 1)
function Base.codeunit(go::GappedOlig, i::Integer)
    @boundscheck 1 <= i <= length(go) || throw(BoundsError(go, i))
    parent_len = length(parent(go))
    cum_gaps = 0
    ungapped_pos = 0 
    gapped_pos = 0
    
    for (start, len) in go.gaps
        positions_before_gap = start - ungapped_pos
        gapped_end = gapped_pos + positions_before_gap
        if gapped_pos < i <= gapped_end
            return UInt8(parent(go)[ungapped_pos + (i - gapped_pos)])
        end
        gapped_pos = gapped_end
        if gapped_pos < i <= gapped_pos + len
            return UInt8('-')
        end
        gapped_pos += len
        cum_gaps += len
        ungapped_pos = start
    end
    
    if gapped_pos < i <= gapped_pos + (parent_len - ungapped_pos)
        return UInt8(parent(go)[ungapped_pos + (i - gapped_pos)])
    end
    throw(BoundsError(go, i))
end

Base.isvalid(::AbstractOlig) = true
Base.isvalid(olig::AbstractOlig, i::Int) = 1 <= i <= length(olig)
Base.length(olig::Olig) = length(String(olig))
Base.length(olig::DegenOlig) = length(String(olig))
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
hasgaps(ov::OligView{DegenOlig}) = false
hasgaps(go::GappedOlig) = !isempty(go.gaps)
hasgaps(ov::OligView) = any(c == '-' for c in ov)

const EMPTY_GAPS_VECTOR = Vector{Pair{Int, Int}}()
getgaps(::AbstractOlig) = EMPTY_GAPS_VECTOR
getgaps(go::GappedOlig) = go.gaps

n_unique_oligs(::AbstractOlig) = BigInt(1)
n_unique_oligs(d::DegenOlig) = d.n_unique_oligs
n_unique_oligs(ov::OligView) = reduce(*, (IUPAC_COUNTS[base] for base in ov), init=BigInt(1))
n_unique_oligs(go::GappedOlig) = n_unique_oligs(parent(go))

n_deg_pos(::AbstractOlig) = 0
n_deg_pos(d::DegenOlig) = d.n_deg_pos
n_deg_pos(ov::OligView) = count(char -> char in DEGEN_BASES, ov)
n_deg_pos(go::GappedOlig) = n_deg_pos(parent(go))

function Base.convert(::Type{T}, o::AbstractOlig) where T<:AbstractOlig
    if T === typeof(o)
        return o
    elseif T <: AbstractGapped
        return GappedOlig(String(o), description(o))
    elseif T <: AbstractDegen
        return DegenOlig(String(o), description(o))
    else T <: AbstractOlig
        return Olig(String(o), description(o))
    end
end

_base_olig_type(::Type{T}) where {T<:Olig} = Olig
_base_olig_type(::Type{T}) where {T<:DegenOlig} = DegenOlig
_base_olig_type(::Type{GappedOlig}) = GappedOlig
_base_olig_type(::Type{OligView{U}}) where {U} = _base_olig_type(U)
_base_olig_type(::Type{T}) where {T<:AbstractOlig} = T  # fallback

_is_degenerate_type(::Type{T}) where {T<:AbstractOlig} = false
_is_degenerate_type(::Type{T}) where {T<:AbstractDegen} = true
_is_degenerate_type(::Type{OligView{T}}) where {T} = _is_degenerate_type(T)

_has_gaps_type(::Type{T}) where {T<:AbstractOlig} = T <: Union{GappedOlig, OligView{<:GappedOlig}}

function Base.promote_rule(::Type{T1}, ::Type{T2}) where {T1<:AbstractOlig, T2<:AbstractOlig}
    is_deg_type = _is_degenerate_type(T1) || _is_degenerate_type(T2)
    is_gapped_type = _has_gaps_type(T1) || _has_gaps_type(T2)
    return is_gapped_type ? GappedOlig : is_deg_type ? DegenOlig : Olig
end
Base.promote_rule(::Type{<:Union{String, SubString}}, ::Type{T}) where {T<:AbstractOlig} = T
Base.promote_rule(::Type{T}, ::Type{<:Union{String, SubString}}) where {T<:AbstractOlig} = T

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
    d = description(parent(iter))
    descr = isempty(d) ? "Non-degen sample" : "Non-degen sample from: $d"
    return Olig(String(buffer), descr), (indices, options, lens, buffer, n)
end
Base.length(iter::NonDegenIterator) = iter.n_variants
Base.eltype(::Type{<:NonDegenIterator}) = Olig

nondegens(olig::Olig) = isempty(olig) ? Tuple{}() : (olig,)
nondegens(go::GappedOlig) = hasgaps(go) ?
    error("Cannot iterate over sequence with gaps") :
    nondegens(DegenOlig(go))

nondegens(deg::DegenOlig) = n_deg_pos(deg) == 0 ?
    nondegens(Olig(deg)) :
    NonDegenIterator(deg, n_unique_oligs(deg))
nondegens(ov::OligView{T}) where T = nondegens(T(ov))

function sampleChar(olig::AbstractOlig)
    n = length(olig)
    n == 0 && throw(ArgumentError("Cannot sample character from empty oligomer"))
    return olig[rand(1:n)]
end

function sampleView(olig::AbstractOlig, len::Int)
    n = length(olig)
    len <= 0 && throw(ArgumentError("Length must be positive, got $len"))
    len > n && throw(ArgumentError("Requested view length $len exceeds oligomer length $n"))
    start = rand(1:(n - len + 1))
    return olig[start:start+len-1]
end

sampleNondeg(o::Olig) = o
sampleNondeg(o::OligView{Olig}) = Olig(o)
function sampleNondeg(d::T) where T <: AbstractOlig
    OutType = _base_olig_type(T)
    isempty(d) && return OutType()
    buffer = Vector{Char}(undef, length(d))
    @inbounds for (i, c) in enumerate(d)
        options = get(IUPAC_B2V, c, ('-',))
        buffer[i] = rand(options)
    end
    d = description(d)
    descr = isempty(d) ? "Non-degen sample" : "Non-degen sample of $d"
    return OutType(String(buffer), descr)
end

function sample_max_gc(d::T) where T <: AbstractOlig
    (isempty(d) || T == Olig) && return d
    
    buffer = Vector{Char}(undef, length(d))
    @inbounds for (i, c) in enumerate(d)
        options = get(MAX_GC_OPTIONS, c, ('-',))
        buffer[i] = rand(options)
    end
    
    base_descr = description(d)
    descr = isempty(base_descr) ? "Max GC content sample" : "Max GC content sample of $base_descr"
    
    return T(String(buffer), descr)
end

function sample_min_gc(d::T) where T <: AbstractOlig
    (isempty(d) || T == Olig) && return d
    
    buffer = Vector{Char}(undef, length(d))
    @inbounds for (i, c) in enumerate(d)
        options = get(MIN_GC_OPTIONS, c, ('-',))
        buffer[i] = rand(options)
    end
    
    base_descr = description(d)
    descr = isempty(base_descr) ? "Min GC content sample" : "Min GC content sample of $base_descr"
    
    return T(String(buffer), descr)
end


Base.iterate(go::GappedOlig) = length(go) == 0 ? nothing : iterate(go, (1, 1, 1, 0))
function Base.iterate(go::GappedOlig, state::NTuple{4, Int})
    pos, ungapped_pos, gap_idx, remaining = state
    if pos > length(go)
        return nothing
    end
    gaps = getgaps(go)
    if remaining > 0
        char = '-'
        new_remaining = remaining - 1
        new_ungapped = ungapped_pos
        new_gap_idx = new_remaining == 0 ? gap_idx + 1 : gap_idx
    elseif gap_idx <= length(gaps) && ungapped_pos == gaps[gap_idx].first
        _, len = gaps[gap_idx]
        char = '-'
        new_remaining = len - 1
        new_ungapped = ungapped_pos
        new_gap_idx = new_remaining == 0 ? gap_idx + 1 : gap_idx
    else
        if ungapped_pos > length(parent(go))
            throw(BoundsError(parent(go), ungapped_pos))
        end
        char = parent(go)[ungapped_pos]
        new_ungapped = ungapped_pos + 1
        new_remaining = 0
        new_gap_idx = gap_idx
    end
    return char, (pos + 1, new_ungapped, new_gap_idx, new_remaining)
end
Base.iterate(olig::AbstractOlig) = length(olig) == 0 ? nothing : (olig[1], 2)
Base.iterate(olig::AbstractOlig, i::Int) = i>length(olig) ? nothing : (olig[i],i+1)


function unfolded_proportion(olig::AbstractOlig; temp::Real=37.0, max_samples::Int=1000)::Float64
    _ext_unfolded_prop(olig; temp=temp, max_samples=max_samples)
end

_ext_unfolded_prop(olig; temp, max_samples) = error(
    "`unfolded_proportion` function requires SeqFold library to be loaded.\n" *
    "In order to get this functionality, please `]add SeqFold` to your project\n" *
    "and load it with `using SeqFold` before constructing primers."
)


include("show_oligs.jl")
end # module
