__precompile__(false)
module SeqFoldExt

using Statistics

using SeqFold

using AlignedOligs.Oligs
using AlignedOligs.Primers


SeqFold.revcomp(olig::T) where T <: AbstractOlig = T(SeqFold.revcomp(String(olig); table=Oligs.DNA_COMP_TABLE_DEG))
SeqFold.complement(olig::T) where T <: AbstractOlig = T(SeqFold.complement(String(olig); table=Oligs.DNA_COMP_TABLE_DEG))

function SeqFold.gc_content(olig::AbstractOlig)::Float64
    if hasgaps(olig)
        return SeqFold.gc_content(parent(olig))  # Ignore gaps
    end
    if isempty(olig)
        return NaN
    end
    total_gc = 0.0
    for c in olig
        total_gc += Oligs.IUPAC_GC_CONTENT[c]
    end
    
    return total_gc / length(olig)
end
function SeqFold.dg(olig::AbstractOlig; temp::Real=37.0, max_variants::Int=1000)::Float64
    hasgaps(olig) && error("Folding not supported for gapped sequences")
    isempty(olig) && return NaN
    if n_unique_oligs(olig) == 1
        return SeqFold.dg(String(olig); temp=temp)
    end
    dgs = Vector{Float64}(undef, min(n_unique_oligs(olig), max_variants))
    @inbounds if n_unique_oligs(olig) > max_variants
        for k in 1:max_variants
            o = sampleNondeg(olig)
            dgs[k] = SeqFold.dg(String(o); temp=temp)
        end
    else
        i = 1
        for o in nondegens(olig)
            dgs[i] = SeqFold.dg(String(o); temp=temp)
            i += 1
        end
    end
    return minimum(dgs)
end
function SeqFold.dg_cache(olig::AbstractOlig; temp::Real = 37.0)::Matrix{Float64}
    if hasgaps(olig)
        error("Free energy cache not supported for gapped sequences")
    else
        return SeqFold.dg_cache(String(Olig(olig)); temp=temp)
    end
end
function SeqFold.tm(
    olig1::AbstractOlig,
    olig2::AbstractOlig; 
    
    conditions=:pcr,
    conf_int::Real=0.8,
    max_variants::Int=1000,
    kwargs...
)
    (hasgaps(olig1) || hasgaps(olig2)) && error("Melting temperature calculation not supported for gapped sequences")
    !(0 < conf_int <= 1) && throw(ArgumentError("conf_int must be in range (0, 1]"))

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
            o1 = sampleNondeg(olig1)
            o2 = sampleNondeg(olig2)
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
SeqFold.tm(
    olig::AbstractOlig;
    
    conditions=:pcr,
    conf_int::Real=0.9,
    max_variants::Int=1000,
    kwargs...
) = SeqFold.tm(olig, SeqFold.complement(olig); conditions=conditions, conf_int=conf_int, max_variants=max_variants, kwargs...)

function SeqFold.tm_cache(olig1::AbstractOlig, olig2::AbstractOlig; conditions=:pcr, kwargs...)::Matrix{Float64}
    (hasgaps(olig1) || hasgaps(olig2)) && error("Melting temperature cache not supported for gapped sequences")
    SeqFold.tm_cache(String(Olig(olig1)), String(Olig(olig2)); conditions=conditions, kwargs...)
end
SeqFold.tm_cache(
    olig::AbstractOlig;
    
    conditions=:pcr,
    kwargs...
)::Matrix{Float64} = SeqFold.tm_cache(olig, SeqFold.complement(olig); conditions=conditions, kwargs...)

SeqFold.dot_bracket(
    olig::AbstractOlig, structs::Vector{SeqFold.Structure}
) = SeqFold.dot_bracket(String(Olig(olig)), structs)

SeqFold.gc_cache(olig::AbstractOlig)::Matrix{Float64} = SeqFold.gc_cache(String(Olig(olig)))


SeqFold.tm(primer::AbstractPrimer) = primer.tm
SeqFold.dg(primer::AbstractPrimer) = primer.dg
SeqFold.gc_content(primer::AbstractPrimer) = primer.gc

Primers._ext_revcomp(args...; kwargs...) = SeqFold.revcomp(args...; kwargs...)
Primers._ext_tm(args...; kwargs...) = SeqFold.tm(args...; kwargs...)
Primers._ext_dg(args...; kwargs...) = SeqFold.dg(args...; kwargs...)
Primers._ext_gc_content(args...; kwargs...) = SeqFold.gc_content(args...; kwargs...)

end  # module
