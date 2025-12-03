module Primers

export AbstractPrimer
export Primer
export construct_primers, best_pairs

using ..Oligs
using ..Alignments

using ProgressMeter
using Statistics

abstract type AbstractPrimer{T<:Union{Olig,DegenOlig}} end

struct Primer{T} <: AbstractPrimer{T}
    msa::AbstractMSA
    pos::UnitRange{Int}
    is_forward::Bool
    consensus::T
    tail_length::Int
    tm::@NamedTuple{mean::Float64, conf::Tuple{Float64, Float64}}
    dg::Float64
    gc::Float64
    slack::Float64
end

function Primer(
    msa::AbstractMSA,
    interval::UnitRange{Int};
    is_forward::Bool=true,
    tail_length::Int=3,
    max_samples::Int=1000,
    tm_conf_int::Real=0.8,
    tm_conds=:pcr,
    dg_temp::Real=37.0,
    slack::Real=0.0,
    descr="Primer for $(nseqs(msa)) seq MSA at positions $interval"
)
    _cons = consensus_degen(msa, interval; slack=slack)
    gapped_cons = is_forward ? _cons : _ext_revcomp(_cons)
    
    underlying_olig = DegenOlig(String(gapped_cons), string(descr))
    
    Tm = _ext_tm(underlying_olig; max_variants=max_samples, conf_int=tm_conf_int, conditions=tm_conds)
    dG = _ext_dg(underlying_olig; max_variants=max_samples, temp=dg_temp)
    GC = _ext_gc_content(underlying_olig)
    Primer(msa, interval, is_forward, underlying_olig, tail_length, Tm, dG, GC, slack)
end

# These are overloaded in ext/SeqFoldExt.jl to load SeqFold.jl library dynamically
_ext_revcomp(args...; kwargs...) = error(
    "Primer construction requires SeqFold library to be loaded.\n" *
    "In order to get this functionality, please `]add SeqFold` to your project\n" *
    "and load it with `using SeqFold` before constructing primers.")
_ext_tm(args...; kwargs...) = error(
    "Primer construction requires SeqFold library to be loaded.\n" *
    "In order to get this functionality, please `]add SeqFold` to your project\n" *
    "and load it with `using SeqFold` before constructing primers.")
_ext_dg(args...; kwargs...) = error(
    "Primer construction requires SeqFold library to be loaded.\n" *
    "In order to get this functionality, please `]add SeqFold` to your project\n" *
    "and load it with `using SeqFold` before constructing primers.")
_ext_gc_content(args...; kwargs...) = error(
    "Primer construction requires SeqFold library to be loaded.\n" *
    "In order to get this functionality, please `]add SeqFold` to your project\n" *
    "and load it with `using SeqFold` before constructing primers."
)

Base.String(primer::AbstractPrimer) = String(primer.consensus)
Base.length(primer::AbstractPrimer) = length(primer.consensus)
Base.isempty(primer::AbstractPrimer) = isempty(primer.consensus)
Base.iterate(primer::AbstractPrimer, state...) = iterate(primer.consensus, state...)
Base.getindex(primer::AbstractPrimer, i::Int) = getindex(primer.consensus, i)
Base.getindex(primer::AbstractPrimer, r::UnitRange{Int}) = getindex(primer.consensus, r)

Base.convert(::Type{DegenOlig}, primer::AbstractPrimer) = primer.consensus

Oligs.n_unique_oligs(primer::AbstractPrimer) = n_unique_oligs(primer.consensus)
Oligs.n_deg_pos(primer::AbstractPrimer) = n_deg_pos(primer.consensus)
Oligs.description(primer::AbstractPrimer) = description(primer.consensus)
Oligs.hasgaps(::AbstractPrimer) = false
Oligs.nondegens(primer::AbstractPrimer) = nondegens(primer.consensus)
Oligs.olig_range(primer::AbstractPrimer) = primer.pos

function construct_primers(
    msa::AbstractMSA;
    is_forward::Bool=true,
    length_range::UnitRange{Int}=17:23,
    tail_length::Int=3,
    tail_degen_pos::Int=0,
    head_degen_pos::Int=5,
    slack::Real=0.05,
    gc_range::UnitRange{Int}=40:60,
    tm_range::UnitRange{Int}=55:60,
    min_delta_g::Real=-5.0,
    min_msadepth::Float64=0.75,
    max_olig_variants::Int=100,
    max_samples::Int=5000,
    tm_conf_int::Real=0.2,
    tm_conds=:pcr,
    dg_temp::Real=mean(tm_range)
)::Vector{Primer}
    0 ≤ slack < 1 || throw(ArgumentError("slack must be in [0,1)"))
    
    primers = Primer[]
    L = length(msa)
    base_count = get_base_count(msa)
    prog = Progress(length(length_range); desc="Constructing... ", color=:white, barlen=10)
    l = ReentrantLock()
    
    Threads.@threads for len in length_range
        len > L && continue
        tail_len = min(tail_length, len)
        head_len = len - tail_len
        head_len < 0 && continue

        for startpos in 1:(L - len + 1)
            interval::UnitRange{Int} = startpos:(startpos + len - 1)
            depths = msadepth(msa, interval)
            any(<(min_msadepth), depths) && continue
            
            if is_forward
                head_interval = startpos:(startpos + head_len - 1)
                tail_interval = (startpos + head_len):(startpos + len - 1)
            else
                head_interval = (startpos + tail_len):(startpos + len - 1)
                tail_interval = startpos:(startpos + tail_len - 1)
            end
            
            if head_len > 0
                head_freqs = @view base_count[:, head_interval]
                head_deg = sum(count(>(slack), col) > 1 for col in eachcol(head_freqs))
                head_deg > head_degen_pos && continue
            end
            
            if tail_len > 0
                tail_freqs = @view base_count[:, tail_interval]
                tail_deg = sum(count(>(slack), col) > 1 for col in eachcol(tail_freqs))
                tail_deg > tail_degen_pos && continue
            end

            _cons = consensus_degen(msa, interval; slack=slack)
            gapped_cons = is_forward ? _cons : _ext_revcomp(_cons)
            hasgaps(gapped_cons) && continue
            
            cons = DegenOlig(gapped_cons)
            n_unique_oligs(cons) > max_olig_variants && continue

            cons = n_unique_oligs(cons) == 1 ? Olig(cons) : cons
            
            gc = _ext_gc_content(cons)
            !(gc_range.start / 100 <= gc <= gc_range.stop / 100) && continue
            
            dg_val = _ext_dg(cons; max_variants=max_samples, temp=dg_temp)
            dg_val < min_delta_g && continue
            
            Tm = _ext_tm(cons; max_variants=max_samples, conf_int=tm_conf_int, conditions=tm_conds)
            (tm_range.stop < first(Tm.conf) || last(Tm.conf) < tm_range.start) && continue
            
            primer = Primer(msa, interval, is_forward, cons, tail_len, Tm, dg_val, gc, slack)

            lock(l)
            try
                push!(primers, primer)
                next!(prog)
            finally
                unlock(l)
            end
        end
    end
    
    return primers
end

function best_pairs(
    forwards::Vector{Primer},
    reverses::Vector{Primer};
    amplicon_len::UnitRange{Int}=0:9999,
    max_tm_diff::Real=4.0
)::Vector{Pair{Primer}}
    pairs = Pair{Primer}[]

    isempty(forwards) || isempty(reverses) && return pairs

    all(p -> p.is_forward, forwards)  || throw(ArgumentError("All forwards must be forward primers"))
    all(p -> !p.is_forward, reverses) || throw(ArgumentError("All reverses must be reverse primers"))
    
    anymsa = root(first(forwards).msa)
    all(root(p.msa) == anymsa for p in forwards) || throw(ArgumentError("All primers must refer to the same MSA"))
    all(root(p.msa) == anymsa for p in reverses) || throw(ArgumentError("All primers must refer to the same MSA"))
    
    @showprogress desc="Matching primer pairs..." enabled=(length(forwards)>1000) barlen=10 for f in forwards
        for r in reverses
            if f.pos.stop >= r.pos.start
                # overlapping primers
                continue
            end
            if max_tm_diff < abs(f.tm.mean - r.tm.mean)
                # too big melting T difference
                continue
            end
                
            amplicon = r.pos.stop - f.pos.start + 1
            if amplicon in amplicon_len
                push!(pairs, f => r)
            end
        end
    end
    
    sort!(pairs; by = p -> abs(p.first.tm.mean - p.second.tm.mean))
    return pairs
end


include("show_primers.jl")

end # module