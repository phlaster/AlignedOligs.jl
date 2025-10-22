struct Primer <: AbstractOlig
    msa::AbstractMSA
    pos::UnitRange{Int}
    is_forward::Bool
    consensus::DegenerateOlig
    tm::@NamedTuple{mean::Float64, conf::Tuple{Float64, Float64}}
    dg::Float64
    gc::Float64

    function Primer(
        msa::AbstractMSA,
        interval::UnitRange{Int};
        is_forward::Bool=true,
        ignore_minors::Real=0.0,
        max_samples::Int=1000,
        tm_conf_int::Real=0.8,
        tm_conditions::Symbol=:pcr,
        dg_temp::Real=37.0
    )
        0 ≤ ignore_minors ≤ 1 || throw(ArgumentError("ignore_minors must be in [0,1]"))
        checkbounds(msa, :, interval)
        gapped_cons = consensus_degen(msa, interval)
        if !is_forward
            gapped_cons = revcomp(gapped_cons)
        end
        dir_str = is_forward ? "forward" : "reverse"
        degeneracy = n_deg_pos(gapped_cons) > 0 ? "Degenerate" : "Non degenerate"
        new_desc = "$degeneracy $dir_str primer for $(nseqs(msa)) seq MSA at positions $interval"
        if ignore_minors > 0
            new_desc *= ", ignore_minors=$ignore_minors"
        end
        underlying_olig = try 
            DegenerateOlig(String(gapped_cons), new_desc)
        catch e
            if e isa ErrorException && occursin("'-'", e.msg)
                throw(ArgumentError("Selected range has all gaps in some positions; cannot construct primer"))
            else
                rethrow()
            end
        end
        Tm = tm(underlying_olig; max_variants=max_samples, conf_int=tm_conf_int, conditions=tm_conditions)
        dG = dg(underlying_olig; max_variants=max_samples, temp=dg_temp)
        GC = SeqFold.gc_content(underlying_olig)
        new(msa, interval, is_forward, underlying_olig, Tm, dG, GC)
    end
end

# String interface delegation
Base.String(primer::Primer) = String(primer.consensus)
Base.length(primer::Primer) = length(primer.consensus)
Base.lastindex(primer::Primer) = lastindex(primer.consensus)
Base.getindex(primer::Primer, i::Int) = primer.consensus[i]
Base.getindex(primer::Primer, r::UnitRange{Int}) = primer.consensus[r]
Base.iterate(primer::Primer) = iterate(primer.consensus)
Base.iterate(primer::Primer, state) = iterate(primer.consensus, state)
Base.ncodeunits(primer::Primer) = ncodeunits(primer.consensus)
Base.codeunit(primer::Primer, i::Integer) = codeunit(primer.consensus, i)
Base.isvalid(primer::Primer, i::Integer) = isvalid(primer.consensus, i)
Base.isempty(primer::Primer) = isempty(primer.consensus)

n_unique_oligs(primer::Primer) = n_unique_oligs(primer.consensus)
n_deg_pos(primer::Primer) = n_deg_pos(primer.consensus)
description(primer::Primer) = description(primer.consensus)
hasgaps(primer::Primer) = false
nondegens(primer::Primer) = nondegens(primer.consensus)
Base.rand(rng::AbstractRNG, primer::Primer) = rand(rng, primer.consensus)
Base.rand(primer::Primer) = rand(primer.consensus)

Base.:(==)(p1::Primer, p2::Primer) = String(p1) == String(p2)
Base.:(==)(p::Primer, s::AbstractString) = String(p) == s
Base.:(==)(s::AbstractString, p::Primer) = p == s


function construct_primers(
    msa::AbstractMSA;
    ignore_minors::Float64=0.0,
    is_forward::Bool=true,
    len_range::UnitRange{Int}=18:22,
    gc_range::UnitRange{Int}=40:60,
    tm_range::UnitRange{Int}=55:60,
    min_delta_g::Float64=-5.0,
    tail_len::Int=3,
    max_degen_body::Int=5,
    max_degen_tail::Int=0,
    max_variants::Int=100,
    max_gap_prop::Float64=0.25
)
    L = length(msa)
    primers = Primer[]
    for l in len_range
        l < tail_len && continue  # tail_len must be <= l
        for start in 1:(L - l + 1)
            range = start:(start + l - 1)
            # Compute average gap proportion
            gap_fracs = [1 - sum(get_base_count(msa, j)) for j in range]
            avg_gap = mean(gap_fracs)
            if avg_gap > max_gap_prop
                continue
            end
            # Build primer
            try
                primer = Primer(msa, range, is_forward, ignore_minors)
                # Check variants
                n_vars = n_unique_oligs(primer)
                if n_vars > max_variants
                    continue
                end
                # Check degenerate positions
                degen_pos = [c in DEGEN_BASES for c in String(primer)]
                body_count = count(degen_pos[1:(l - tail_len)])
                tail_count = count(degen_pos[(l - tail_len + 1):end])
                if body_count > max_degen_body || tail_count > max_degen_tail
                    continue
                end
                # Compute stats
                stats = olig_stats(primer)
                gc_perc = round(Int, stats.gc * 100)
                if !(gc_range.start <= gc_perc <= gc_range.stop)
                    continue
                end
                if !(tm_range.start <= stats.Tm <= tm_range.stop)
                    continue
                end
                if stats.dg < min_delta_g
                    continue
                end
                push!(primers, primer)
            catch e
                if e isa ArgumentError && occursin("all gaps", e.msg)
                    continue
                else
                    rethrow()
                end
            end
        end
    end
    return primers
end