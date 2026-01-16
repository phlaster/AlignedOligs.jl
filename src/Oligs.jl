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

    function Olig(seq::AbstractString, descr::Union{AbstractString,Integer})
        if !isempty(seq) && !isvalid(Olig, seq)
            invalid_chars = join(setdiff(Set(seq), NON_DEGEN_BASES), ", ")
            error("Olig contains unallowed characters: $invalid_chars")
        end
        
        seq = uppercase(seq)
        descr_str = string(descr)
        
        new(seq, descr_str)
    end
end
struct DegenOlig <: AbstractDegen
    seq::String
    n_deg_pos::Int
    n_unique_oligs::BigInt
    description::String

    function DegenOlig(seq::AbstractString, descr::Union{AbstractString,Integer})
        if !isempty(seq) && !isvalid(DegenOlig, seq)
            invalid_chars = join(setdiff(Set(seq), ALL_BASES), ", ")
            error("DegenOlig contains unallowed characters: $invalid_chars")
        end
        
        seq = uppercase(seq)
        descr_str = string(descr)

        n_deg = count(char -> char in DEGEN_BASES, seq)
        n_unique = reduce(*, (IUPAC_COUNTS[char] for char in seq), init=BigInt(1))
        
        new(seq, n_deg, n_unique, descr_str)
    end
end
struct GappedOlig <: AbstractGapped
    parent::DegenOlig
    gaps::Vector{Pair{Int, Int}}
    total_length::Int

    function GappedOlig(seq::AbstractString, descr::Union{AbstractString,Integer})
        if !isempty(seq) && !isvalid(GappedOlig, seq)
            invalid_chars = join(setdiff(Set(seq), BASES_W_GAPS), ", ")
            error("DegenOlig contains unallowed characters: $invalid_chars")
        end
        
        seq = uppercase(seq)
        total_len = length(seq)
        descr_str = string(descr)
        parent_seq = filter(!=('-'), seq)
        parent_olig = DegenOlig(parent_seq, descr_str)
        
        gaps = Pair{Int,Int}[]
        parent_pos = 0
        i = 1
        
        @inbounds while i <= total_len
            if seq[i] != '-'
                parent_pos += 1
                i += 1
            else
                gap_start = parent_pos + 1
                len = 0
                while i <= total_len && seq[i] == '-'
                    len += 1
                    i += 1
                end
                push!(gaps, gap_start => len)
            end
        end

        return new(parent_olig, gaps, total_len)
    end
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
    NonDegenIterator(o::T) where T<:AbstractOlig = new{T}(o, n_unique_oligs(o))
end


##########################
#      Base methods      #
##########################
#### NonDegenIterator ####
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
Base.parent(ndi::NonDegenIterator) = ndi.olig




######### Oligs ##########
Base.:(==)(o::AbstractOlig, s::SubString{<:Base.AnnotatedString}) = String(o) == String(s)
Base.:(==)(s::SubString{<:Base.AnnotatedString}, o::AbstractOlig) = o == s
Base.:(==)(o::AbstractOlig, s::Base.AnnotatedString) = String(o) == String(s)
Base.:(==)(s::Base.AnnotatedString, o::AbstractOlig) = o == s
Base.:(==)(a::GappedOlig, b::GappedOlig) = (parent(a) == parent(b)) && (getgaps(a) == getgaps(b)) && (length(a) == length(b))

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

Base.uppercase(olig::AbstractOlig) = String(olig)
Base.String(olig::AbstractOlig) = olig.seq
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

Base.parent(olig::AbstractOlig) = olig
Base.parent(olig::OligView) = olig.parent
Base.parent(olig::GappedOlig) = olig.parent

Base.ncodeunits(olig::AbstractOlig) = length(olig)

Base.codeunit(::AbstractOlig) = UInt8
Base.codeunit(olig::AbstractOlig, i::Integer) = codeunit(String(olig), i)
Base.codeunit(ov::OligView, i::Integer) = codeunit(parent(ov), olig_range(ov).start + i - 1)
function Base.codeunit(go::GappedOlig, i::Integer)
    @boundscheck 1 <= i <= length(go) || throw(BoundsError(go, i))
    parent_len = length(parent(go))
    ungapped_pos = 1
    gapped_pos = 0
    
    for (start, len) in go.gaps
        num_parent_chars_before_gap = (start - 1) - (ungapped_pos - 1)
        
        if gapped_pos < i <= gapped_pos + num_parent_chars_before_gap
            offset = i - gapped_pos
            return UInt8(parent(go)[ungapped_pos + offset - 1])
        end
        
        gapped_pos += num_parent_chars_before_gap
        ungapped_pos += num_parent_chars_before_gap
        
        if gapped_pos < i <= gapped_pos + len
            return UInt8('-')
        end
        gapped_pos += len
    end
    
    offset = i - gapped_pos
    if ungapped_pos + offset - 1 > parent_len
        throw(BoundsError(go, i))
    end
    return UInt8(parent(go)[ungapped_pos + offset - 1])
end

Base.isvalid(::AbstractOlig) = true
Base.isvalid(olig::AbstractOlig, i::Int) = 1 <= i <= length(olig)
Base.isvalid(::Type{Olig}, s::AbstractString) = all(c -> uppercase(c) in NON_DEGEN_BASES, s)
Base.isvalid(::Type{DegenOlig}, s::AbstractString) = all(c -> uppercase(c) in ALL_BASES, s)
Base.isvalid(::Type{GappedOlig}, s::AbstractString) = all(c -> uppercase(c) in BASES_W_GAPS, s)

Base.length(olig::AbstractOlig) = length(String(olig))
Base.length(go::GappedOlig) = go.total_length
Base.length(ov::OligView) = length(olig_range(ov))
Base.isempty(olig::AbstractOlig) = length(olig) == 0
Base.lastindex(olig::AbstractOlig) = length(olig)

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


const EMPTY_OLIG = Olig("", "")
Olig() = EMPTY_OLIG
Olig(olig::Olig) = olig
Olig(seq::AbstractString) = Olig(seq, description(seq))

const EMPTY_DEGENERATE = DegenOlig("", "")
DegenOlig() = EMPTY_DEGENERATE
DegenOlig(olig::DegenOlig) = olig
DegenOlig(seq::AbstractString) = DegenOlig(seq, description(seq))

const EMPTY_GAPPED = GappedOlig(DegenOlig(), "")
GappedOlig() = EMPTY_GAPPED
GappedOlig(olig::GappedOlig) = olig
GappedOlig(seq::AbstractString) = GappedOlig(seq, description(seq))



##########################
#      Type methods      #
##########################

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
getgaps(ov::OligView{GappedOlig}) = ov |> GappedOlig |> getgaps

n_unique_oligs(::AbstractOlig) = BigInt(1)
n_unique_oligs(d::DegenOlig) = d.n_unique_oligs
n_unique_oligs(ov::OligView) = reduce(*, (IUPAC_COUNTS[base] for base in ov), init=BigInt(1))
n_unique_oligs(go::GappedOlig) = n_unique_oligs(parent(go))

n_deg_pos(::AbstractOlig) = 0
n_deg_pos(d::DegenOlig) = d.n_deg_pos
n_deg_pos(ov::OligView) = count(char -> char in DEGEN_BASES, ov)
n_deg_pos(go::GappedOlig) = n_deg_pos(parent(go))

nondegens(olig::Olig) = isempty(olig) ? Tuple{}() : (olig,)
nondegens(go::GappedOlig) = hasgaps(go) ?
    error("Cannot iterate over sequence with gaps") :
    nondegens(DegenOlig(go))

nondegens(deg::DegenOlig) = n_deg_pos(deg) == 0 ?
    nondegens(Olig(deg)) :
    NonDegenIterator(deg)
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


_base_olig_type(::Type{T}) where {T<:Olig} = Olig
_base_olig_type(::Type{T}) where {T<:DegenOlig} = DegenOlig
_base_olig_type(::Type{GappedOlig}) = GappedOlig
_base_olig_type(::Type{OligView{U}}) where {U} = _base_olig_type(U)
_base_olig_type(::Type{T}) where {T<:AbstractOlig} = T  # fallback

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

function unfolded_proportion(args...; kwargs...)
    error(
        "`unfolded_proportion` function requires SeqFold library to be loaded.\n" *
        "In order to get this functionality, please `]add SeqFold` to your project\n" *
        "and load it with `using SeqFold`."
    )
end


include("show_oligs.jl")
end # module
