abstract type AbstractMSA end

struct MSA <: AbstractMSA
    seqs::Vector{GappedOlig}
    base_count::Matrix{Float64}  # 4 x L, proportions of A=1, C=2, G=3, T=4

        function MSA(seqs::Vector{GappedOlig})
        if isempty(seqs)
            return new(seqs, zeros(4, 0))
        end

        L = length(first(seqs))
        all(length(s) == L for s in seqs) || throw(ArgumentError("All gapped sequences must have the same length"))

        n = length(seqs)
        base_count = zeros(4, L)

        for j in 1:L
            counts = zeros(Int, 4)
            num_non_gaps = 0
            for s in seqs
                c = s[j]
                if c != '-'
                    idx = findfirst(==(c), "ACGT")
                    if !isnothing(idx)
                        counts[idx] += 1
                        num_non_gaps += 1
                    end
                end
            end
            if num_non_gaps > 0
                base_count[:, j] = counts / num_non_gaps
            end
        end

        return new(seqs, base_count)
    end
end

struct MSAView <: AbstractMSA
    parent::AbstractMSA
    rows::UnitRange{Int}  # absolute indices relative to root MSA
    cols::UnitRange{Int}  # absolute indices relative to root MSA
end

function root(msa::MSA)
    return msa
end

function root(v::MSAView)
    return root(v.parent)
end

function MSA(fasta::String; aligned::Bool=true, gap_tolerance::Float64=0.5)
    0.0 < gap_tolerance < 1.0 || throw(ArgumentError("gap_tolerance must be in (0,1)"))

    seqs_desc = []
    FastaReader(fasta) do fr
        for (i, (desc, seq)) in enumerate(fr)
            desc = isempty(desc) ? "seq$i" : desc
            push!(seqs_desc, (seq, desc))
        end
    end

    if !aligned
        mktemp() do input_path, io
            writefasta(io, seqs_desc; check_description=false)
            aligned_output = read(`$(MAFFT_jll.mafft()) $input_path`, String)
            seqs_desc = []
            # Read directly from IOBuffer
            FastaReader(IOBuffer(aligned_output)) do fr
                for (desc, seq) in fr
                    push!(seqs_desc, (seq, desc))
                end
            end
        end
    end

    gapped_oligs = [GappedOlig(seq, desc) for (seq, desc) in seqs_desc]
    msa = MSA(gapped_oligs)

    # Apply gap_tolerance to filter low-abundance nucleotides
    if gap_tolerance > 0
        L = width(msa)
        new_seqs = Vector{GappedOlig}(undef, nseqs(msa))
        for i in 1:nseqs(msa)
            seq_chars = Char[]
            for j in 1:L
                p = msa.base_count[:, j]
                c = msa.seqs[i][j]
                if c != '-' && sum(p) > 0
                    idx = findfirst(==('A'), "ACGT")  # Example, replace with actual base logic
                    if p[idx] < gap_tolerance
                        push!(seq_chars, '-')  # Replace low-abundance bases with gaps
                    else
                        push!(seq_chars, c)
                    end
                else
                    push!(seq_chars, c)
                end
            end
            new_seqs[i] = GappedOlig(String(seq_chars), description(msa.seqs[i]))
        end
        return MSA(new_seqs)
    end

    return msa
end

function MSA(fasta::String; aligned::Bool=true)
    seqs_desc = Tuple{String,String}[]
    FastaReader(fasta) do fr
        for (i, (desc, seq)) in enumerate(fr)
            desc = isempty(desc) ? "seq$i" : desc
            push!(seqs_desc, (desc, seq))
        end
    end

    if !aligned
        output_path = fasta * ".aln"
        mktemp() do input_path, _
            open(input_path, "w") do io
                FastaWriter(io) do fw
                    for (desc, seq) in seqs_desc
                        writeentry(fw, desc, seq)
                    end
                end
            end
            aligned_output = read(`$(MAFFT_jll.mafft()) $input_path`, String)
            open(output_path, "w") do io
                write(io, aligned_output)
            end
        end
        seqs_desc = Tuple{String,String}[]
        FastaReader(output_path) do fr
            for (desc, seq) in fr
                push!(seqs_desc, (seq, desc))
            end
        end
    end

    # Validate sequences
    for (desc, seq) in seqs_desc
        allowed = aligned ? "ACGT-" : "ACGT"
        upper_seq = uppercase(seq)
        if !all(c -> c in allowed, upper_seq)
            invalid_chars = setdiff(unique(upper_seq), collect(allowed))
            throw(ArgumentError("Sequence '$desc' contains invalid characters: $(join(invalid_chars, ", ")). Only $(join(collect(allowed), ", ")) allowed."))
        end
    end

    gapped_oligs = [GappedOlig(seq, desc) for (seq, desc) in seqs_desc]
    return MSA(gapped_oligs)
end

# Slicing
function Base.getindex(msa::AbstractMSA, rows::UnitRange{Int}, cols::UnitRange{Int})
    root_msa = root(msa)
    abs_rows = if isa(msa, MSA)
        rows
    else
        msa.rows.first + rows.first - 1 : msa.rows.first + rows.last - 1
    end
    abs_cols = if isa(msa, MSA)
        cols
    else
        msa.cols.first + cols.first - 1 : msa.cols.first + cols.last - 1
    end
    return MSAView(root_msa, abs_rows, abs_cols)
end

Base.getindex(msa::AbstractMSA, rows::Colon, cols::UnitRange{Int}) = msa[1:nseqs(msa), cols]
Base.getindex(msa::AbstractMSA, rows::UnitRange{Int}, cols::Colon) = msa[rows, 1:width(msa)]
Base.getindex(msa::AbstractMSA, rows::Colon, cols::Colon) = msa[1:nseqs(msa), 1:width(msa)]

# Accessors
nseqs(msa::MSA) = length(msa.seqs)
nseqs(v::MSAView) = length(v.rows)

width(msa::MSA) = size(msa.base_count, 2)
width(v::MSAView) = length(v.cols)

function getsequence(msa::MSA, row::Int)
    return msa.seqs[row]
end

function getsequence(v::MSAView, row::Int)
    abs_row = v.rows.first + row - 1
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
    abs_pos = v.cols.first + pos - 1
    return get_base_count(root(v), abs_pos)
end

function get_base_count(msa::AbstractMSA)
    L = width(msa)
    counts = zeros(4, L)
    for j in 1:L
        counts[:, j] = get_base_count(msa, j)
    end
    return counts
end

# Consensus functions
function consensus_major(msa::AbstractMSA, pos::Int)
    p = get_base_count(msa, pos)
    if sum(p) == 0
        return '-'
    end
    return "ACGT"[argmax(p)]
end

function consensus_degen(msa::AbstractMSA, pos::Int)
    p = get_base_count(msa, pos)
    if sum(p) == 0
        return '-'
    end
    active = findall(>(0), p)
    if isempty(active)
        return '-'
    end
    bs = sort(collect("ACGT"[active]))
    return get(IUPAC_V2B, bs, 'N')  # fallback to 'N' if not found
end

# Full consensus sequences
function consensus_major(msa::AbstractMSA)
    return join(consensus_major(msa, j) for j in 1:width(msa))
end

function consensus_degen(msa::AbstractMSA)
    return join(consensus_degen(msa, j) for j in 1:width(msa))
end

# Placeholder for dry_msa (remove all-gap columns)
function dry_msa(msa::AbstractMSA)
    non_gap_cols = [j for j in 1:width(msa) if sum(get_base_count(msa, j)) > 0]
    if non_gap_cols isa UnitRange
        return msa[:, non_gap_cols]
    else
        # Materialize if non-contiguous
        seqs = [GappedOlig(join(getsequence(msa, i)[non_gap_cols])) for i in 1:nseqs(msa)]
        return MSA(seqs)
    end
end
