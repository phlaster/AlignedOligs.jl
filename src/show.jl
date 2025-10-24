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
    max_display_height = max(5, terminal_height - 5)
    max_display_width = max(20, terminal_width - 10)

    if seq_length <= max_display_width
        println(io, ":")
        if n_sequences <= max_display_height
            for i in 1:n_sequences
                println(io, getsequence(msa, i))
            end
        else
            n_first = max_display_height ÷ 2
            n_last = max_display_height - n_first - 1
            n_first = min(n_first, n_sequences)
            n_last = min(n_last, n_sequences)
            n_last = min(n_last, n_sequences - n_first)

            for i in 1:n_first
                println(io, getsequence(msa, i))
            end
            if n_first + n_last < n_sequences
                println(io, "⋅⋅⋅")
            end
            start_idx_last = max(n_first + 1, n_sequences - n_last + 1)
            for i in start_idx_last:n_sequences
                if i > 0 && i <= n_sequences
                    println(io, getsequence(msa, i))
                end
            end
        end
    else
        println(io, ":")
        half_width = max_display_width ÷ 2 - 2
        half_width = max(1, half_width)

        if n_sequences <= max_display_height
            for i in 1:n_sequences
                seq_str = string(getsequence(msa, i))
                seq_len = length(seq_str)
                if seq_len <= max_display_width
                    println(io, seq_str)
                else
                    start_end_idx = max(1, min(half_width, seq_len))
                    end_start_idx = max(start_end_idx + 1, seq_len - half_width + 1)
                    start_part = seq_str[1:start_end_idx]
                    end_part = seq_str[end_start_idx:end]
                    println(io, start_part * "⋅⋅⋅" * end_part)
                end
            end
        else
            n_first = max_display_height ÷ 2
            n_last = max_display_height - n_first - 1
            n_first = min(n_first, n_sequences)
            n_last = min(n_last, n_sequences)
            n_last = min(n_last, n_sequences - n_first)

            for i in 1:n_first
                 seq_str = string(getsequence(msa, i))
                 seq_len = length(seq_str)
                 if seq_len <= max_display_width
                     println(io, seq_str)
                 else
                     start_end_idx = max(1, min(half_width, seq_len))
                     end_start_idx = max(start_end_idx + 1, seq_len - half_width + 1)
                     start_part = seq_str[1:start_end_idx]
                     end_part = seq_str[end_start_idx:end]
                     println(io, start_part * "⋅⋅⋅" * end_part)
                 end
            end

            if n_first + n_last < n_sequences
                println(io, " " ^ half_width * " ⋮ ")
            end

            start_idx_last = max(n_first + 1, n_sequences - n_last + 1)
            for i in start_idx_last:n_sequences
                seq_str = string(getsequence(msa, i))
                seq_len = length(seq_str)
                if seq_len <= max_display_width
                    println(io, seq_str)
                else
                    start_end_idx = max(1, min(half_width, seq_len))
                    end_start_idx = max(start_end_idx + 1, seq_len - half_width + 1)
                    start_part = seq_str[1:start_end_idx]
                    end_part = seq_str[end_start_idx:end]
                    println(io, start_part * "⋅⋅⋅" * end_part)
                end
            end
        end
    end
end