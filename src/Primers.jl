struct Primer <: AbstractOlig
    msa::AbstractMSA
    pos::UnitRange{Int}
    is_forward::Bool
    consensus::DegenerateOlig
    tail_length::Int
    tm::@NamedTuple{mean::Float64, conf::Tuple{Float64, Float64}}
    dg::Float64
    gc::Float64
end

function Primer(
    msa::AbstractMSA,
    interval::UnitRange{Int};
    is_forward::Bool=true,
    tail_length::Int=3,
    max_samples::Int=1000,
    tm_conf_int::Real=0.8,
    tm_conds::Symbol=:pcr,
    dg_temp::Real=37.0
)
    thresh = min_tresh(msa)
    gapped_cons = consensus_degen(msa, interval)
    if !is_forward
        gapped_cons = revcomp(gapped_cons)
    end
    dir_str = is_forward ? "forward" : "reverse"
    degeneracy = n_deg_pos(gapped_cons) > 0 ? "Degenerate" : "Non degenerate"
    new_desc = "$degeneracy $dir_str primer for $(nseqs(msa)) seq MSA at positions $interval"
    if thresh > 0
        new_desc *= ", minor_thresh=$thresh"
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
    Tm = tm(underlying_olig; max_variants=max_samples, conf_int=tm_conf_int, conditions=tm_conds)
    dG = dg(underlying_olig; max_variants=max_samples, temp=dg_temp)
    GC = SeqFold.gc_content(underlying_olig)
    Primer(msa, interval, is_forward, underlying_olig, tail_length, Tm, dG, GC)
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
Base.rand(rng::AbstractRNG, primer::Primer) = rand(rng, primer.consensus)
Base.rand(primer::Primer) = rand(primer.consensus)
Base.:(==)(p1::Primer, p2::Primer) = String(p1) == String(p2)
Base.:(==)(p::Primer, s::AbstractString) = String(p) == s
Base.:(==)(s::AbstractString, p::Primer) = p == s

n_unique_oligs(primer::Primer) = n_unique_oligs(primer.consensus)
n_deg_pos(primer::Primer) = n_deg_pos(primer.consensus)
description(primer::Primer) = description(primer.consensus)
hasgaps(primer::Primer) = false
nondegens(primer::Primer) = nondegens(primer.consensus)
olig_range(primer::Primer) = primer.pos

function construct_primers(
    msa::AbstractMSA;
    is_forward::Bool=true,
    length_range::UnitRange{Int}=17:23,
    tail_length::Int=3,
    tail_degen_pos::Int=0,
    head_degen_pos::Int=5,
    head_slack::Float64=0.15,
    tail_slack::Float64=0.05,
    gc_range::UnitRange{Int}=40:60,
    tm_range::UnitRange{Int}=55:60,
    min_delta_g::Real=-5.0,
    min_msadepth::Float64=0.75,
    max_olig_variants::Int=100,
    max_samples::Int=5000,
    tm_conf_int::Real=0.2,
    tm_conds::Symbol=:pcr,
    dg_temp::Real=mean(tm_range)
)::Vector{Primer}
    primers = Primer[]
    L = length(msa)
    base_count = get_base_count(msa)
    prog = Progress(length(length_range); desc="Constructing... ", color=:white, barlen=10)
    Threads.@threads for len in length_range
        len > L && continue
        tail_len = min(tail_length, len)
        head_len = len - tail_len
        head_len < 0 && continue
        for startpos in 1:(L - len + 1)
            rng = startpos:(startpos + len - 1)
            depths = msadepth(msa, rng)
            any(<(min_msadepth), depths) && continue
            if is_forward
                head_rng = startpos:(startpos + head_len - 1)
                tail_rng = (startpos + head_len):(startpos + len - 1)
            else
                head_rng = (startpos + tail_len):(startpos + len - 1)
                tail_rng = startpos:(startpos + tail_len - 1)
            end
            if head_len > 0
                head_freqs = @view base_count[:, head_rng]
                head_deg = sum(count(>(head_slack), col) > 1 for col in eachcol(head_freqs))
                head_deg > head_degen_pos && continue
            end
            if tail_len > 0
                tail_freqs = @view base_count[:, tail_rng]
                tail_deg = sum(count(>(tail_slack), col) > 1 for col in eachcol(tail_freqs))
                tail_deg > tail_degen_pos && continue
            end
            gapped_cons = consensus_degen(msa, rng)
            if !is_forward
                gapped_cons = SeqFold.revcomp(gapped_cons)
            end
            underlying_olig = try
                DegenerateOlig(String(gapped_cons), "Primer for $(height(msa))seq MSA")
            catch e
                if isa(e, ErrorException) && occursin("'-'", e.msg)
                    continue
                else
                    rethrow()
                end
            end
            # APPLY FILTERS
            n_unique_oligs(underlying_olig) > max_olig_variants && continue
            
            gc = SeqFold.gc_content(underlying_olig)
            !(gc_range.start / 100 <= gc <= gc_range.stop / 100) && continue
            
            dg_val = SeqFold.dg(underlying_olig; max_variants=max_samples, temp=dg_temp)
            dg_val < min_delta_g && continue
            
            Tm = SeqFold.tm(underlying_olig; max_variants=max_samples, conf_int=tm_conf_int, conditions=tm_conds)
            (tm_range.stop < first(Tm.conf) || last(Tm.conf) < tm_range.start) && continue
            
            primer = Primer(msa, rng, is_forward, underlying_olig, tail_len, Tm, dg_val, gc)

            l = ReentrantLock()
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

SeqFold.tm(primer::Primer) = primer.tm
SeqFold.dg(primer::Primer) = primer.dg
SeqFold.gc_content(primer::Primer) = primer.gc

function best_pairs(
    forwards::Vector{Primer},
    reverses::Vector{Primer};
    amplicon_len::UnitRange{Int}=0:9999
)::Vector{Pair{Primer}}
    pairs = Pair{Primer}[]

    isempty(forwards) || isempty(reverses) && return pairs

    all(p -> p.is_forward, forwards)  || throw(ArgumentError("All forwards must be forward primers"))
    all(p -> !p.is_forward, reverses) || throw(ArgumentError("All reverses must be reverse primers"))
    anymsa = root(rand(forwards).msa)
    all(root(p.msa) == anymsa for p in forwards) || throw(ArgumentError("All primers must refer to the same MSA"))
    all(root(p.msa) == anymsa for p in reverses) || throw(ArgumentError("All primers must refer to the same MSA"))
    
    @showprogress desc="Matching primer pairs..." barlen=10 for f in forwards
        for r in reverses
            if f.pos.stop >= r.pos.start
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