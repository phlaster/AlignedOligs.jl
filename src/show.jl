function _truncate_seq(seq::AbstractString, max_width::Int=20)
    if length(seq) > max_width
        return seq[1:max(0, max_width-3)] * "..."
    else
        return seq
    end
end

function Base.show(io::IO, olig::Olig)
    seq_display = _truncate_seq(String(olig))
    print(io, "Olig(\"", seq_display, "\", len=", length(olig))
    if !isempty(description(olig))
        print(io, ", desc=\"", description(olig), "\"")
    end
    print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", olig::Olig)
    println(io, "Olig")
    println(io, "  Sequence: ", String(olig))
    println(io, "  Length: ", length(olig))
    print(io, "  Description: ")
    if isempty(description(olig))
        print(io, "(none)")
    else
        print(io, "\"", description(olig), "\"")
    end
end

function Base.show(io::IO, deg::DegenerateOlig)
    seq_display = _truncate_seq(String(deg))
    print(io, "DegenerateOlig(\"", seq_display, "\", len=", length(deg))
    print(io, ", n_deg=", n_deg_pos(deg))
    vars_str = n_unique_oligs(deg) > 10000 ? ">10k" : string(n_unique_oligs(deg))
    print(io, ", vars=", vars_str)
    if !isempty(description(deg))
        print(io, ", desc=\"", description(deg), "\"")
    end
    print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", deg::DegenerateOlig)
    println(io, "DegenerateOlig")
    println(io, "  Sequence: ", String(deg))
    println(io, "  Length: ", length(deg))
    println(io, "  Degenerate positions: ", n_deg_pos(deg))
    println(io, "  Unique variants: ", n_unique_oligs(deg))
    print(io, "  Description: ")
    if isempty(description(deg))
        print(io, "(none)")
    else
        print(io, "\"", description(deg), "\"")
    end
end

function Base.show(io::IO, go::GappedOlig)
    seq_display = _truncate_seq(String(go))
    print(io, "GappedOlig(\"", seq_display, "\", len=", length(go))
    print(io, ", gaps=", length(go.gaps))
    if !isempty(description(go))
        print(io, ", desc=\"", description(go), "\"")
    end
    print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", go::GappedOlig)
    println(io, typeof(go))
    println(io, "  Gapped sequence: ", String(go))
    println(io, "  Length (with gaps): ", length(go))
    println(io, "  Underlying Olig: ", parent(go))
    println(io, "  Gaps: ", length(go.gaps))
    print(io, "  Description: ")
    if isempty(description(go))
        print(io, "(none)")
    else
        print(io, "\"", description(go), "\"")
    end
end

function Base.show(io::IO, ov::OligView)
    seq_display = _truncate_seq(String(ov))
    print(io, "OligView(\"", seq_display, "\", len=", length(ov))
    print(io, ", range=", ov.range)
    if !isempty(description(ov))
        print(io, ", desc=\"", description(ov), "\"")
    end
    print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", ov::OligView)
    println(io, "OligView{", typeof(parent(ov)), "}")
    println(io, "  Viewed sequence: ", String(ov))
    println(io, "  Length: ", length(ov))
    println(io, "  Range: ", ov.range)
    print(io, "  Parent description: ", description(ov))
end

function Base.show(io::IO, primer::Primer)
    max_width = 15
    seq = String(primer)
    if length(seq) > max_width
        seq_display = seq[1:max_width-3] * "..."
    else
        seq_display = seq
    end
    dir_str = primer.is_forward ? "forward" : "reverse"
    print(io, "Primer(\"", seq_display, "\", len=$(length(primer)), pos=$(primer.pos), $dir_str")
    print(io, ", degen=$(n_deg_pos(primer)), variants=$(n_unique_oligs(primer) > 10000 ? ">10k" : n_unique_oligs(primer))")
    print(io, ", Tm=$(round(primer.tm.mean, digits=1))°C, ΔG=$(round(primer.dg, digits=2))kcal/mol, GC=$(round(primer.gc * 100, digits=1))%)")
end

function Base.show(io::IO, ::MIME"text/plain", primer::Primer)
    dir_str = primer.is_forward ? "Forward" : "Reverse"
    ndeg = n_deg_pos(primer)
    deg_status = ndeg > 0 ?
        "degenerate primer with $ndeg deg. positions" :
        "non-degenerate primer"
    println(io, "$dir_str $deg_status")

    L = length(primer.msa)
    pos = primer.pos
    s, e = pos.start, pos.stop
    term_width = displaysize(io)[2] - 1
    L_str = string(L)
    bar_start = "|"
    bar_end = L_str * "|"
    eq_len = max(0, term_width - length(bar_start) - length(bar_end))
    bar = bar_start * repeat('=', eq_len) * bar_end

    if L == 1
        map_start_col = length(bar_start) + 1
        map_end_col = map_start_col + eq_len - 1
        col_s = col_e = map_start_col
    else
        map_start_col = length(bar_start) + 1
        map_end_col = map_start_col + eq_len - 1
        scale = (eq_len - 1) / (L - 1.0)
        col_s = map_start_col + round(Int, (s - 1) * scale)
        col_e = map_start_col + round(Int, (e - 1) * scale)
    end

    if primer.is_forward
        three_prime_col = col_e
        label = "\\" * string(s) * ":" * string(e) * ">"
    else
        three_prime_col = col_s
        label = "<" * string(s) * ":" * string(e) * "\\"
    end

    label_len = length(label)
    indent = max(0, three_prime_col - label_len)
    label_line = " "^indent * label

    if primer.is_forward
        println(io, label_line)
        println(io, bar)
        println(io)
    else
        println(io)
        println(io, bar)
        println(io, label_line)
    end

    println(io, "  Sequence: ", String(primer))
    println(io, "  Length: ", length(primer))
    println(io, "  Positions: ", olig_range(primer))
    println(io, "  Unique variants: ", n_unique_oligs(primer))
    println(io, "  Melting temperature: ", round(primer.tm.mean, digits=1), "°C (", round(primer.tm.conf[1], digits=1), "⋅", round(primer.tm.conf[2], digits=1), "°C)")
    println(io, "  Min ΔG: ", round(primer.dg, digits=2), " kcal/mol")
    println(io, "  GC content: ", round(primer.gc * 100, digits=1), "%")
    print(io, "  Description: \"", description(primer), "\"")
end

function Base.show(io::IO, ::MIME"text/plain", pp::Pair{Primer})
    fwd, rev = pp.first, pp.second
    
    if !fwd.is_forward || rev.is_forward
        invoke(show, Tuple{IO, MIME"text/plain", Pair}, io, MIME"text/plain"(), pp)
        return
    end
    
    try
        msa = root(fwd.msa)
        if msa !== root(rev.msa)
            invoke(show, Tuple{IO, MIME"text/plain", Pair}, io, MIME"text/plain"(), pp)
            return
        end
        
        N = nseqs(msa)
        L = length(msa)
        amp_start = fwd.pos.start
        amp_end = rev.pos.stop
        amp_len = amp_end - amp_start + 1
        overlap = fwd.pos.stop >= rev.pos.start ? "!!! OVERLAPPING !!! " : ""
        header = "$(overlap)PCR primer pair for $N seq. MSA, amplicon: $amp_start:$amp_end ($(amp_len)bp)"
        println(io, header)

        term_width = displaysize(io)[2] - 1
        bar_start = "|"
        bar_end = string(L) * "|"
        eq_len = max(0, term_width - length(bar_start) - length(bar_end))
        bar = bar_start * repeat('=', eq_len) * bar_end
        scale = L == 1 ? 0.0 : (eq_len - 1) / (L - 1.0)
        map_start_col = length(bar_start) + 1
        col_amp_start = map_start_col + (L > 1 ? round(Int, (amp_start - 1) * scale) : 0)
        col_amp_end = map_start_col + (L > 1 ? round(Int, (amp_end - 1) * scale) : 0)
        arrow_line = [' ' for _ in 1:term_width]
        if col_amp_start <= term_width && col_amp_start >= 1
            arrow_line[col_amp_start] = '>'
        end
        if col_amp_end <= term_width && col_amp_end >= 1 && col_amp_end != col_amp_start
            arrow_line[col_amp_end] = '<'
        end
        for i in (col_amp_start + 1):(col_amp_end - 1)
            if 1 <= i <= term_width
                arrow_line[i] = '_'
            end
        end
        label = fwd.pos.stop >= rev.pos.start ? "" : "$(amp_len)bp"
        label_len = length(label)
        inner_len = col_amp_end - col_amp_start - 1
        if inner_len >= label_len + 2
            mid_col = col_amp_start + div(col_amp_end - col_amp_start + 1, 2)
            label_start = mid_col - div(label_len, 2)
            if label_start >= col_amp_start + 2 && label_start + label_len - 1 <= col_amp_end - 2
                for (j, c) in enumerate(label)
                    pos = label_start + j - 1
                    if 1 <= pos <= term_width
                        arrow_line[pos] = c
                    end
                end
            end
        end
        println(io, join(arrow_line))
        println(io, bar)
        println(io, "Forward: ", String(fwd), " at ", fwd.pos.start, ":", fwd.pos.stop)
        println(io, "Reverse: ", String(rev), " at ", rev.pos.start, ":", rev.pos.stop)
        mean_tm = (fwd.tm.mean + rev.tm.mean) / 2
        delta_tm = abs(fwd.tm.mean - rev.tm.mean) / 2
        print(io, "Tm: ", round(mean_tm, digits=1), "±", round(delta_tm, digits=1), " °C")
    catch e
        invoke(show, Tuple{IO, MIME"text/plain", Pair}, io, MIME"text/plain"(), pp)
    end
end

const BASE_COLORS = Dict{Char, Symbol}(
    'A' => :green,
    'C' => :blue,
    'G' => :yellow,
    'T' => :red,
    'M' => :light_cyan,   # A/C
    'R' => :light_green,  # A/G
    'W' => :light_green,  # A/T
    'S' => :light_blue,   # C/G
    'Y' => :light_magenta,# C/T
    'K' => :light_red,    # G/T
    'V' => :cyan,         # A/C/G
    'H' => :magenta,      # A/C/T
    'D' => :green,        # A/G/T
    'B' => :blue,         # C/G/T
    'N' => :white,
    '-' => :normal
)

function Base.show(io::IO, msa::AbstractMSA)
    n_sequences = nseqs(msa)
    if n_sequences == 0
        print(io, "Empty ", typeof(msa))
        return
    end
    seq_length = length(msa)
    print(io, typeof(msa), " with $n_sequences sequences of length $seq_length")
    if seq_length == 0
        return
    end

    local terminal_height, terminal_width
    try
        terminal_height, terminal_width = displaysize(io)
    catch
        terminal_height, terminal_width = 24, 80
    end
    max_display_height = max(5, terminal_height - 7)
    n_display_seqs = min(n_sequences, max_display_height)

    # Collect and process descriptions for displayed sequences
    processed_descs = Vector{String}(undef, n_display_seqs)
    has_desc = false
    for i in 1:n_display_seqs
        desc = description(getsequence(msa, i))
        processed = replace(replace(desc, '\n' => ' '), '\t' => ' ')
        processed_descs[i] = processed
        if !isempty(processed)
            has_desc = true
        end
    end

    desc_width = 0
    padded_descs = String[]
    if has_desc
        max_possible_desc = min(50, floor(Int, 0.2 * terminal_width))
        max_desc_len = maximum(length, processed_descs)
        desc_width = min(max_possible_desc, max_desc_len)
        padded_descs = Vector{String}(undef, n_display_seqs)
        for i in 1:n_display_seqs
            trunc = processed_descs[i][1:min(end, desc_width)]
            padded_descs[i] = rpad(trunc, desc_width)
        end
        desc_width += 2  # for separator
    end

    max_display_width = max(20, terminal_width - desc_width - 7)

    needs_width_ellipsis = seq_length > max_display_width
    max_seq_chars = needs_width_ellipsis ? max_display_width - 3 : max_display_width
    displayed_cols_range = 1:min(seq_length, max_seq_chars)

    # Compute absolute columns if MSAView
    abs_cols = if msa isa MSAView
        msa.cols.start .+ (displayed_cols_range .- 1)
    else
        displayed_cols_range
    end

    # Compute bold_vec for displayed columns
    bold_vec = Vector{Bool}(undef, length(displayed_cols_range))
    if _is_full_height(msa)
        for dj in 1:length(displayed_cols_range)
            det = msadet(msa, dj)
            bold_vec[dj] = det < 1.0 && det > 0.0
        end
    else
        for dj in 1:length(displayed_cols_range)
            pos_counts = zeros(Float64, 4)
            for i in 1:n_sequences
                c = getsequence(msa, i, dj)
                probs = IUPAC_PROBS[c]
                pos_counts .+= probs
            end
            s = sum(pos_counts)
            bold_vec[dj] = false
            if s > 0.0
                det = maximum(pos_counts) / s
                bold_vec[dj] = det < 1.0
            end
        end
    end

    println(io, ":")

    for i in 1:n_display_seqs
        if has_desc
            print(io, padded_descs[i], " >")
        end
        for dj in 1:length(displayed_cols_range)
            c = getsequence(msa, i, dj)
            col = get(BASE_COLORS, c, :normal)
            bold = bold_vec[dj]
            printstyled(io, c; color=col, bold=bold, reverse=true)
        end
        if needs_width_ellipsis
            printstyled(io, "..."; color=:light_black, reverse=true)
        end
        println(io)
    end

    if n_display_seqs < n_sequences
        println(io, "...")
    end
end