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

struct OligView{T<:AbstractOlig} <: AbstractOlig
    parent::T
    range::UnitRange{Int}

    function OligView{T}(parent, range) where T<:AbstractOlig
        new(parent(parent), range)
    end
end

Base.String(olig::Olig) = olig.seq
Base.String(olig::DegenerateOlig) = olig.seq
Base.String(ov::OligView) = String(ov.parent)[ov.range]

parent(olig::AbstractOlig) = olig
parent(olig::OligView) = olig.parent

Base.ncodeunits(olig::AbstractOlig) = ncodeunits(String(olig))
Base.codeunit(olig::AbstractOlig) = codeunit(String(olig))
Base.codeunit(olig::AbstractOlig, i::Integer) = codeunit(String(olig), i)
Base.isvalid(olig::AbstractOlig, i::Integer) = 1 <= i <= length(olig) && isvalid(String(olig), i)
Base.iterate(olig::AbstractOlig) = iterate(String(olig))
Base.iterate(olig::AbstractOlig, state::Integer) = iterate(String(olig), state)
Base.lastindex(olig::AbstractOlig) = length(olig)
Base.length(olig::AbstractOlig) = length(String(olig))
Base.isempty(olig::AbstractOlig) = isempty(String(olig))

description(olig::Olig) = olig.description
description(olig::DegenerateOlig) = olig.description
description(ov::OligView) = description(ov.parent)

Olig(chars::Vector{Char}, description::AbstractString = "") = Olig(String(chars), description)

function Olig(olig::AbstractOlig)
    if n_deg_pos(olig) > 0
        throw(InexactError(:Olig, Olig, olig))
    end
    return Olig(String(olig), description(olig))
end

Olig() = EMPTY_OLIG

function DegenerateOlig(seq::AbstractString, description = "")
    isempty(seq) && return DegenerateOlig()

    seq = uppercase(seq)
    seq_chars = Set(seq)
    if !issubset(seq_chars, ALL_BASES)
        error("DegenerateOlig contains unallowed characters: $(join(setdiff(seq_chars, ALL_BASES), ", "))")
    end
    
    n_degenerate = count(char -> char in DEGEN_BASES, seq)
    n_possible = reduce(*, IUPAC_COUNTS[char] for char in seq, init=BigInt(1))
    return DegenerateOlig(seq, n_degenerate, n_possible, string(description))
end

DegenerateOlig() = EMPTY_OLIG

DegenerateOlig(chars::Vector{Char}, description::AbstractString = "") = DegenerateOlig(String(chars), description)
DegenerateOlig(olig::AbstractOlig) = DegenerateOlig(String(olig), n_deg_pos(olig), n_unique_oligs(olig), description(olig))

function Base.getindex(olig::T, r::UnitRange{Int}) where {T<:AbstractOlig}
    @boundscheck checkbounds(1:lastindex(olig), r)
    OligView(olig, r)
end

n_unique_oligs(::Olig) = 1
n_unique_oligs(d::DegenerateOlig) = d.n_unique_oligs
n_unique_oligs(::OligView{Olig}) = 1
n_unique_oligs(ov::OligView{DegenerateOlig}) = reduce(*, IUPAC_COUNTS[String(ov)[i]] for i in 1:length(ov), init=BigInt(1))

n_deg_pos(::Olig) = 0
n_deg_pos(d::DegenerateOlig) = d.n_deg_pos
n_deg_pos(::OligView{Olig}) = 0
n_deg_pos(ov::OligView{DegenerateOlig}) = count(char -> char in DEGEN_BASES, String(ov))

function Base.:*(o1::AbstractOlig, o2::AbstractOlig)
    if n_deg_pos(o1) > 0 || n_deg_pos(o2) > 0
        return DegenerateOlig(String(o1) * String(o2), n_deg_pos(o1)+n_deg_pos(o2), n_unique_oligs(o1)*n_unique_oligs(o2), "concatenated")
    else
        return Olig(String(o1) * String(o2), "concatenated")
    end
end

Base.convert(::Olig, o::AbstractOlig) = Olig(o)
Base.convert(::DegenerateOlig, o::AbstractOlig) = DegenerateOlig(o)
Base.promote_rule(::Type{Olig}, ::Type{DegenerateOlig}) = DegenerateOlig

function Base.show(io::IO, olig::AbstractOlig)
    color = get(io, :color, false)
    yellow = color ? "\e[33m" : ""
    reset = color ? "\e[0m" : ""

    seq = String(olig)
    max_width = 20
    if length(seq) > max_width
        seq_display = seq[1:min(max_width-3, end)] * "..."
        seq_info = " $yellow$(length(seq)) nt$reset"
    else
        seq_display = seq
        seq_info = ""
    end
    
    if olig isa DegenerateOlig
        deg_str = n_deg_pos(olig) == 0 ? "non-degen" :
                  n_deg_pos(olig) == 1 ? "single-degen" :
                  "$(n_deg_pos(olig)) deg"
        
        variants = n_unique_oligs(olig)
        variants_str = variants > 10_000 ? ">10k" : string(variants)
        print(io, "DegenerateOlig: $seq_display$seq_info, $deg_str, $variants_str vars")
    elseif olig isa OligView
        print(io, "OligView: $seq_display$seq_info")
    else
        print(io, "Olig: $seq_display$seq_info")
    end
end

function Base.show(io::IO, ::MIME"text/plain", olig::AbstractOlig)
    color = get(io, :color, false)
    yellow = color ? "\e[33m" : ""
    reset = color ? "\e[0m" : ""

    seq = String(olig)
    seq_info = length(seq) > 40 ? " $yellow$(length(seq)) bases total$reset" : ""
    seq_display = length(seq) > 40 ? seq[1:min(40, end)] * "..." : seq
    
    if olig isa DegenerateOlig
        deg_str = n_deg_pos(olig) == 0 ? "non-degenerate" :
                  n_deg_pos(olig) == 1 ? "1 deg" :
                  "$(n_deg_pos(olig)) deg"
        
        variants = n_unique_oligs(olig)
        variants_str = variants > 10_000 ? ">10k" : string(variants)
        
        println(io, "DegenerateOlig: $seq_display$seq_info")
        println(io, "  • $deg_str, $variants_str unique variants")
    elseif olig isa OligView
        println(io, "OligView: $seq_display$seq_info")
    else
        println(io, "Olig: $seq_display$seq_info")
    end
    
    desc = description(olig)
    if !isempty(desc)
        desc_info = length(desc) > 60 ? " $yellow$(length(desc)) bytes total$reset" : ""
        desc_display = length(desc) > 42 ? desc[1:min(42, end)] * "..." : desc
        print(io, "  • $desc_display$desc_info")
    end
end

struct NonDegenIterator{T<:AbstractOlig}
    olig::T
    n_variants::Integer
end

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

    return Olig(String(buffer), description(iter.olig)), (indices, options, lens, buffer, n)
end

Base.length(iter::NonDegenIterator) = iter.n_variants
Base.eltype(::Type{<:NonDegenIterator}) = Olig

nondegens(olig::AbstractOlig) = NonDegenIterator(olig, n_unique_oligs(olig))

function Base.rand(rng::AbstractRNG, olig::Olig)
    return olig  # Non-degenerate, return as-is
end

function Base.rand(rng::AbstractRNG, olig::DegenerateOlig)
    isempty(olig) && return EMPTY_OLIG
    buffer = Vector{Char}(undef, length(olig))
    @inbounds for (i, c) in enumerate(olig.seq)
        options = IUPAC_B2V[c]
        buffer[i] = rand(rng, options)
    end
    return Olig(String(buffer), description(olig))
end

function Base.rand(rng::AbstractRNG, ov::OligView)
    parent_olig = rand(rng, ov.parent)
    return Olig(String(parent_olig)[ov.range], description(ov))
end

Base.rand(olig::AbstractOlig) = rand(Random.GLOBAL_RNG, olig)


# SeqFold specifications
SeqFold.revcomp(olig::T) where T <: AbstractOlig = T(SeqFold.revcomp(String(olig); table=DNA_COMP_TABLE_DEG))
SeqFold.complement(olig::T) where T <: AbstractOlig = T(SeqFold.complement(String(olig); table=DNA_COMP_TABLE_DEG))

function SeqFold.gc_content(olig::AbstractOlig)::Float64
    if isempty(olig)
        return NaN
    end
    total_gc = 0.0
    for c in olig
        total_gc += IUPAC_GC_CONTENT[c]
    end
    
    return total_gc / length(olig)
end

SeqFold.fold(olig::AbstractOlig; temp::Real = 37.0)::Vector{Structure} = SeqFold.fold(String(Olig(olig)); temp=temp)
SeqFold.dg(olig::Olig; temp::Real = 37.0)::Float64 = SeqFold.dg(String(olig); temp=temp)
SeqFold.dg(deg_olig::DegenerateOlig; temp::Real = 37.0)::Float64 = minimum(olig->dg(olig; temp=temp), nondegens(deg_olig))
SeqFold.dg_cache(olig::AbstractOlig; temp::Real = 37.0)::Matrix{Float64} = SeqFold.dg_cache(String(Olig(olig)); temp=temp)
SeqFold.dot_bracket(olig::AbstractOlig, structs::Vector{SeqFold.Structure}) = SeqFold.dot_bracket(String(Olig(olig)), structs)

function SeqFold.tm(olig1::AbstractOlig, olig2::AbstractOlig; conditions=:pcr, conf_int::Real=0.8, max_variants::Int=10000, kwargs...)
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
    SeqFold.tm_cache(String(Olig(olig1)), String(Olig(olig2)); conditions=conditions, kwargs...)
end

function SeqFold.tm_cache(olig::AbstractOlig; conditions=:pcr, kwargs...)::Matrix{Float64}
    SeqFold.tm_cache(olig, SeqFold.complement(olig); conditions=conditions, kwargs...)
end

SeqFold.gc_cache(olig::AbstractOlig)::Matrix{Float64} = SeqFold.gc_cache(String(Olig(olig)))

function Base.stat(olig::AbstractOlig)
    gc = SeqFold.gc_content(olig)
    Tm = tm(olig)
    ΔG = dg(olig)

    return (Tm=Tm, gc=gc, dg=ΔG)
end