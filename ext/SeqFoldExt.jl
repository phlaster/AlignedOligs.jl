module SeqFoldExt

using Statistics

using SeqFold

using AlignedOligs.Oligs
using AlignedOligs.Primers


SeqFold.revcomp(olig::T) where T <: AbstractOlig = T(SeqFold.revcomp(String(olig); table=Oligs.DNA_COMP_TABLE_DEG))
SeqFold.complement(olig::T) where T <: AbstractOlig = T(SeqFold.complement(String(olig); table=Oligs.DNA_COMP_TABLE_DEG))

function SeqFold.gc_content(olig::AbstractOlig)::Float64
    hasgaps(olig) && return SeqFold.gc_content(parent(olig))
    isempty(olig) && return NaN
    total_gc = sum(get(Oligs.IUPAC_GC_CONTENT, c, 0.5) for c in olig, init=0.0)
    return total_gc / length(olig)
end

function logsumexp(x::Vector{T}) where T
    m = maximum(x)
    isinf(m) && return m
    s = sum(exp(val - m) for val in x, init=zero(T))
    return m + log(s)
end

function SeqFold.dg(olig::AbstractOlig; temp::Real=37.0, max_samples::Int=1000, mode::Symbol=:average)::Float64
    hasgaps(olig) && error("Folding not supported for gapped sequences")
    isempty(olig) && return NaN
    total_variants = n_unique_oligs(olig)
    total_variants == 1 && return SeqFold.dg(String(olig); temp=temp)

    T_K::Float64 = temp + 273.15
    R::Float64 = 1.9872e-3
    RT::Float64 = R * T_K
    N::Int = clamp(total_variants, 1, max_samples)
    ΔGs::Float64 = Inf64
    if mode === :average
        log_terms = Vector{Float16}(undef, N)
        if total_variants <= max_samples
            for (i, o) in enumerate(nondegens(olig))
                ΔG = SeqFold.dg(String(o); temp=temp)
                log_terms[i] = -ΔG / RT
            end
        else
            for i in 1:N
                o = sampleNondeg(olig)
                ΔG = SeqFold.dg(String(o); temp=temp)
                log_terms[i] = -ΔG / RT
            end
        end
        log_Q = logsumexp(log_terms) - log(N)
        ΔGs = -RT * log_Q
    elseif mode === :worstcase
        if total_variants <= max_samples
            for o in nondegens(olig)
                ΔGs = min(SeqFold.dg(String(o); temp=temp), ΔGs)
            end
        else
            for _ in 1:max_samples
                o = sampleNondeg(olig)
                ΔGs = min(SeqFold.dg(String(o); temp=temp), ΔGs)
            end
        end
    else
        error("Invalid mode: $mode. Supported modes are :average and :worstcase.")
    end
    return round(ΔGs, digits=2)
end

function SeqFold.dg_cache(olig::AbstractOlig; temp::Real=37.0)::Matrix{Float64}
    hasgaps(olig) && error("Free energy cache not supported for gapped sequences")
    return SeqFold.dg_cache(String(Olig(olig)); temp=temp)
end
function SeqFold.tm(
    olig1::AbstractOlig,
    olig2::AbstractOlig; 
    
    conditions=:pcr,
    conf_int::Real=0.8,
    max_samples::Int=1000,
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
            conf = (mean_Tm, mean_Tm),
            min = mean_Tm,
            max = mean_Tm
        )
    end

    T = zeros(Float64, min(i * j, max_samples))
    Tm_min = -Inf64
    Tm_max = Inf64
    if i * j > max_samples
        for k in 1:max_samples
            o1 = sampleNondeg(olig1)
            o2 = sampleNondeg(olig2)
            @inbounds T[k] = SeqFold.tm(String(o1), String(o2); conditions=conditions, kwargs...)
        end
        Tm_min_sample, Tm_max_sample = extrema(T)
        Tm_min = min(SeqFold.tm(String(sample_min_gc(olig1)), String(sample_min_gc(olig2)); conditions=conditions, kwargs...), Tm_min_sample)
        Tm_max = max(SeqFold.tm(String(sample_max_gc(olig1)), String(sample_max_gc(olig2)); conditions=conditions, kwargs...), Tm_max_sample)
    else
        counter = 1
        for o1 in nondegens(olig1)
            for o2 in nondegens(olig2)
                @inbounds T[counter] = SeqFold.tm(String(o1), String(o2); conditions=conditions, kwargs...)
                counter += 1
            end
        end
        Tm_min, Tm_max = extrema(T)
    end

    mean_Tm = mean(T)
    alpha = (1 - conf_int) / 2
    low = quantile(T, alpha)
    high = quantile(T, 1 - alpha)

    return (
        mean = round(mean_Tm, digits=1),
        conf = (round(low, digits=1), round(high, digits=1)),
        min = Tm_min,
        max = Tm_max
    )
end
SeqFold.tm(
    olig::AbstractOlig;
    
    conditions=:pcr,
    conf_int::Real=0.9,
    max_samples::Int=1000,
    kwargs...
) = SeqFold.tm(olig, SeqFold.complement(olig); conditions=conditions, conf_int=conf_int, max_samples=max_samples, kwargs...)

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

function Oligs.unfolded_proportion(olig; temp, max_samples)
    hasgaps(olig) && error("Folding not supported for gapped sequences")
    isempty(olig) && return NaN
    total_variants = n_unique_oligs(olig)
    T_K::Float64 = temp + 273.15
    R = 1.9872e-3
    RT = R * T_K
    if total_variants == 1
        ΔG = SeqFold.dg(String(olig); temp=temp)
        K_f = exp(-ΔG / RT)
        return clamp(inv(1 + K_f), 0.0, 1.0)
    end
    N = clamp(total_variants, 1, max_samples)
    unfolded_fractions = Vector{Float64}(undef, N)
    if total_variants <= max_samples
        idx = 1
        for o in nondegens(olig)
            ΔG = SeqFold.dg(String(o); temp=temp)
            K_f = exp(-ΔG / RT)
            unfolded_fractions[idx] = inv(1 + K_f)
            idx += 1
        end
        avg_unfolded = mean(unfolded_fractions)
    else
        for i in 1:N
            o = sampleNondeg(olig)
            ΔG = SeqFold.dg(String(o); temp=temp)
            K_f = exp(-ΔG / RT)
            unfolded_fractions[i] = 1 / (1 + K_f)
        end
        avg_unfolded = mean(unfolded_fractions)
    end
    return clamp(avg_unfolded, 0.0, 1.0)
end
function Primers._ext_revcomp(o::AbstractOlig)
    SeqFold.revcomp(o)
end
function Primers._ext_tm(olig::AbstractOlig; max_samples, conf_int, conditions)
    SeqFold.tm(olig; max_samples=max_samples, conf_int=conf_int, conditions=conditions)
end
function Primers._ext_dg(olig::AbstractOlig; max_samples, temp)
    SeqFold.dg(olig; max_samples=max_samples, temp=temp)
end
function Primers._ext_gc_content(olig::AbstractOlig)
    SeqFold.gc_content(olig)
end

end  # module
