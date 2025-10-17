abstract type AbstractMSA end

struct MSA <: AbstractMSA
    seqs::Vector{GappedOlig}
    base_count::Matrix{Float64}
    bootstrap::Int

    function MSA(seqs::Vector{<:GappedOlig}; gap_tolerance::Real=0.0, bootstrap::Int=0, seed::Union{Int,Nothing}=nothing)
        0 ≤ gap_tolerance < 1 || throw(ArgumentError("gap_tolerance must be in [0,1)"))
        bootstrap >= 0 || throw(ArgumentError("bootstrap must be non-negative"))
        isnothing(seed) || Random.seed!(seed)

        if isempty(seqs)
            return new(seqs, zeros(4, 0), bootstrap)
        end

        L = length(first(seqs))
        all(length(s) == L for s in seqs) || throw(ArgumentError("All gapped sequences must have the same length"))
        n = length(seqs)

        if bootstrap > 0
            # Initialize with zeros - only need 4×L memory
            base_count = zeros(4, L)
            
            # For each bootstrap replicate
            for b in 1:bootstrap
                boot_rows = rand(1:n, n)
                boot_seqs = seqs[boot_rows]
                
                # Process one bootstrap replicate
                for j in 1:L
                    pos_counts = zeros(4)
                    num_valid = 0
                    for s in boot_seqs
                        c = s[j]
                        if c != '-'
                            num_valid += 1
                            probs = get(IUPAC_PROBS, c, zeros(4))
                            pos_counts .+= probs
                        end
                    end
                    
                    if num_valid > 0
                        current_freq = pos_counts / num_valid
                        
                        # Online update of the mean
                        base_count[:, j] += (current_freq - base_count[:, j]) / b
                    end
                end
            end
        else
            # Original non-bootstrapped base_count
            base_count = zeros(4, L)
            for j in 1:L
                counts = zeros(4)
                num_valid = 0
                for s in seqs
                    c = s[j]
                    if c != '-'
                        num_valid += 1
                        probs = get(IUPAC_PROBS, c, zeros(4))
                        counts .+= probs
                    end
                end
                if num_valid > 0
                    base_count[:, j] = counts / num_valid
                end
            end
        end

        if gap_tolerance > 0
            new_seqs = Vector{GappedOlig}(undef, n)
            for i in 1:n
                seq_chars = Char[]
                for j in 1:L
                    p = base_count[:, j]
                    c = seqs[i][j]
                    if c != '-' && sum(p) > 0
                        idx = findfirst(==(c), "ACGT")
                        if !isnothing(idx) && p[idx] < gap_tolerance
                            push!(seq_chars, '-')
                        else
                            push!(seq_chars, c)
                        end
                    else
                        push!(seq_chars, c)
                    end
                end
                new_seqs[i] = GappedOlig(String(seq_chars), description(seqs[i]))
            end
            # Recurse with same bootstrap value
            return MSA(new_seqs; gap_tolerance=0.0, bootstrap=bootstrap, seed=seed)
        end

        return new(seqs, base_count, bootstrap)
    end
end

struct MSAView <: AbstractMSA
    parent::AbstractMSA
    rows::UnitRange{Int}
    cols::UnitRange{Int}
end

function root(msa::MSA)
    return msa
end

function root(v::MSAView)
    return root(v.parent)
end

function MSA(fasta::AbstractString; aligned::Bool=true, gap_tolerance::Real=0.0, bootstrap::Int=0)
    seqs_desc = Tuple{String, String}[]
    FastaReader(fasta) do fr
        for (i, (desc, seq)) in enumerate(fr)
            desc = isempty(desc) ? "seq$i" : desc
            push!(seqs_desc, (desc, seq))
        end
    end

    if !aligned
        mktemp() do input_path, io
            FastaWriter(io) do fw
                for (desc, seq) in seqs_desc
                    writeentry(fw, desc, seq)
                end
            end
            close(io)
            aligned_output = read(`$(MAFFT_jll.mafft()) $input_path`, String)
            seqs_desc = Tuple{String, String}[]
            FastaReader(IOBuffer(aligned_output)) do fr
                for (desc, seq) in fr
                    push!(seqs_desc, (desc, seq))
                end
            end
        end
    end

    # Validate sequences
    allowed = aligned ? (NON_DEGEN_BASES..., '-', 'N') : (NON_DEGEN_BASES..., 'N')
    for (desc, seq) in seqs_desc
        upper_seq = uppercase(seq)
        if !all(c -> c in allowed, upper_seq)
            invalid_chars = setdiff(unique(upper_seq), collect(allowed))
            throw(ArgumentError("Sequence '$desc' contains invalid characters: $(join(invalid_chars, ", ")). Only $(join(collect(allowed), ", ")) allowed."))
        end
    end

    gapped_oligs = [GappedOlig(seq, desc) for (desc, seq) in seqs_desc]
    return MSA(gapped_oligs; gap_tolerance=gap_tolerance, bootstrap=bootstrap)
end

function MSA(predicate::Function, fasta::AbstractString; aligned::Bool=true, gap_tolerance::Real=0.0, bootstrap::Int=0)
    seqs_desc = Tuple{String, String}[]
    FastaReader(fasta) do fr
        for (i, (desc, seq)) in enumerate(fr)
            if predicate(desc)
                desc = isempty(desc) ? "seq$i" : desc
                push!(seqs_desc, (seq, desc))
            end
        end
    end

    if !aligned
        mktemp() do input_path, io
            FastaWriter(io) do fw
                for (desc, seq) in seqs_desc
                    writeentry(fw, desc, seq)
                end
            end
            close(io)
            aligned_output = read(`$(MAFFT_jll.mafft()) $input_path`, String)
            seqs_desc = Tuple{String, String}[]
            FastaReader(IOBuffer(aligned_output)) do fr
                for (desc, seq) in fr
                    push!(seqs_desc, (seq, desc))
                end
            end
        end
    end

    # Validate sequences
    allowed = aligned ? (NON_DEGEN_BASES..., '-', 'N') : (NON_DEGEN_BASES..., 'N')
    for (seq, desc) in seqs_desc
        upper_seq = uppercase(seq)
        if !all(c -> c in allowed, upper_seq)
            invalid_chars = setdiff(unique(upper_seq), collect(allowed))
            throw(ArgumentError("Sequence '$desc' contains invalid characters: $(join(invalid_chars, ", ")). Only $(join(collect(allowed), ", ")) allowed."))
        end
    end

    gapped_oligs = [GappedOlig(seq, desc) for (seq, desc) in seqs_desc]
    return MSA(gapped_oligs; gap_tolerance=gap_tolerance, bootstrap=bootstrap)
end

function Base.getindex(msa::AbstractMSA, rows::UnitRange{Int}, cols::UnitRange{Int})
    root_msa = root(msa)
    abs_rows = isa(msa, MSA) ? rows : msa.rows.start + rows.start - 1 : msa.rows.start + rows.stop - 1
    abs_cols = isa(msa, MSA) ? cols : msa.cols.start + cols.start - 1 : msa.cols.start + cols.stop - 1
    return MSAView(root_msa, abs_rows, abs_cols)
end
Base.getindex(msa::AbstractMSA, rows::Colon, cols::UnitRange{Int}) = msa[1:nseqs(msa), cols]
Base.getindex(msa::AbstractMSA, rows::UnitRange{Int}, cols::Colon) = msa[rows, 1:length(msa)]
Base.getindex(msa::AbstractMSA, rows::Colon, cols::Colon) = msa[1:nseqs(msa), 1:length(msa)]

nseqs(msa::MSA) = length(msa.seqs)
nseqs(v::MSAView) = length(v.rows)

Base.length(msa::MSA) = size(msa.base_count, 2)
Base.length(v::MSAView) = length(v.cols)

# aliases
width(msa::AbstractMSA) = length(msa)
height(msa::AbstractMSA) = nseqs(msa)

function getsequence(msa::MSA, row::Int)
    return msa.seqs[row]
end

function getsequence(v::MSAView, row::Int)
    abs_row = v.rows.start + row - 1
    parent_seq = getsequence(root(v), abs_row)
    return parent_seq[v.cols]
end

function getsequence(msa::AbstractMSA, row::Int, col::Int)
    return getsequence(msa, row)[col]
end

function get_base_count(msa::MSA, pos::Int)
    return msa.base_count[:, pos]
end

function get_base_count(v::MSAView, pos::Int)
    abs_pos = v.cols.start + pos - 1
    return get_base_count(root(v), abs_pos)
end

function get_base_count(msa::AbstractMSA)
    L = length(msa)
    counts = zeros(4, L)
    for j in 1:L
        counts[:, j] = get_base_count(msa, j)
    end
    return counts
end

function consensus_major(msa::AbstractMSA, pos::Int; ignore_minors::Real=0.0)
    0 ≤ ignore_minors ≤ 1 || throw(ArgumentError("ignore_minors must be in [0,1]"))
    p = get_base_count(msa, pos)
    if sum(p) == 0
        return '-'
    end
    p_filtered = copy(p)
    for i in 1:4
        if p_filtered[i] ≤ ignore_minors
            p_filtered[i] = 0.0
        end
    end
    if sum(p_filtered) == 0
        return '-'
    end
    return "ACGT"[argmax(p_filtered)]
end

function consensus_degen(msa::AbstractMSA, pos::Int; ignore_minors::Real=0.0)
    0 ≤ ignore_minors ≤ 1 || throw(ArgumentError("ignore_minors must be in [0,1]"))
    p = get_base_count(msa, pos)
    if sum(p) == 0
        return '-'
    end
    active = findall(x -> x > ignore_minors, p)
    if isempty(active)
        return '-'
    end
    bs = sort(collect("ACGT"[active]))
    return get(IUPAC_V2B, bs, 'N')  # fallback to 'N' if not found
end

function consensus_major(msa::AbstractMSA; ignore_minors::Real=0.0)
    seq = join(consensus_major(msa, j; ignore_minors=ignore_minors) for j in 1:length(msa))
    desc = "Major consensus for $(nseqs(msa)) seq MSA"
    if ignore_minors > 0
        desc *= ", ignore_minors=$ignore_minors"
    end
    return GappedOlig(seq, desc)
end

function consensus_degen(msa::AbstractMSA; ignore_minors::Real=0.0)
    seq = join(consensus_degen(msa, j; ignore_minors=ignore_minors) for j in 1:length(msa))
    desc = "Degenerate consensus for $(nseqs(msa)) seq MSA"
    if ignore_minors > 0
        desc *= ", ignore_minors=$ignore_minors"
    end
    return GappedOlig(seq, desc)
end

function dry_msa(msa::AbstractMSA; gap_content::Real=1.0)
    0 ≤ gap_content ≤ 1 || throw(ArgumentError("gap_content must be in [0,1]"))
    non_gap_cols = [j for j in 1:length(msa) if sum(get_base_count(msa, j)) > 0]
    if isempty(non_gap_cols)
        return msa[:, 1:0]
    end
    num_cols = length(non_gap_cols)
    kept_rows = Int[]
    for i in 1:nseqs(msa)
        gap_count = sum((1 for j in non_gap_cols if getsequence(msa, i, j) == '-'), init=0)
        prop = num_cols > 0 ? gap_count / num_cols : 0.0
        if prop < gap_content
            push!(kept_rows, i)
        end
    end
    if isempty(kept_rows)
        return MSA(GappedOlig[]; bootstrap=msa.bootstrap)
    end
    new_seqs = Vector{GappedOlig}(undef, length(kept_rows))
    for (k, row) in enumerate(kept_rows)
        sub_str = join(getsequence(msa, row, j) for j in non_gap_cols)
        new_seqs[k] = GappedOlig(sub_str, description(getsequence(msa, row)))
    end
    return MSA(new_seqs; bootstrap=msa.bootstrap)
end