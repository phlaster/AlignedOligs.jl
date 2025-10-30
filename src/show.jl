function _show_primer_common(io::IO, primer::AbstractPrimer)
    println(io, "  Sequence: ", consensus_sequence(primer))
    println(io, "  Length: ", length(consensus_oligomer(primer)))
    println(io, "  Positions: ", position_range(primer))
    println(io, "  Unique variants: ", n_unique_oligs(primer))
    println(io, "  Melting temperature: ", round(melting_temperature(primer).mean, digits=1), "°C (",
            round(melting_temperature(primer).conf[1], digits=1), "⋅",
            round(melting_temperature(primer).conf[2], digits=1), "°C)")
    println(io, "  Min ΔG: ", round(free_energy(primer), digits=2), " kcal/mol")
    println(io, "  GC content: ", round(gc_content(primer) * 100, digits=1), "%")
    print(io, "  Description: \"", description(primer), "\"")
end

function Base.show(io::IO, primer::AbstractPrimer)
    max_width = 15
    seq = consensus_sequence(primer)
    seq_display = length(seq) > max_width ? seq[1:max_width-3] * "..." : seq

    dir_str = is_forward_primer(primer) ? "forward" : "reverse"
    print(io, "Primer(\"", seq_display, "\", len=$(length(consensus_oligomer(primer))), ",
          "pos=$(position_range(primer).start):$(position_range(primer).stop), $dir_str")

    print(io, ", degen=$(n_deg_pos(primer)), variants=$(n_unique_oligs(primer) > 10000 ? ">10k" : n_unique_oligs(primer)))")
    print(io, ", Tm=$(round(melting_temperature(primer).mean, digits=1))°C, ",
          "ΔG=$(round(free_energy(primer), digits=2))kcal/mol, ",
          "GC=$(round(gc_content(primer) * 100, digits=1))%)")
end

function Base.show(io::IO, ::MIME"text/plain", primer::AbstractPrimer)
    dir_str = is_forward_primer(primer) ? "Forward" : "Reverse"
    ndeg = n_deg_pos(primer)
    deg_status = ndeg > 0 ?
        "degenerate primer with $ndeg deg. positions" :
        "non-degenerate primer"
    println(io, "$dir_str $deg_status")

    L = length(msa_reference(primer))
    pos = position_range(primer)
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

    if is_forward_primer(primer)
        three_prime_col = col_e
        label = "\\" * string(s) * ":" * string(e) * ">"
    else
        three_prime_col = col_s
        label = "<" * string(s) * ":" * string(e) * "\\"
    end

    label_len = length(label)
    indent = max(0, three_prime_col - label_len)
    label_line = " "^indent * label

    if is_forward_primer(primer)
        println(io, label_line)
        println(io, bar)
        println(io)
    else
        println(io)
        println(io, bar)
        println(io, label_line)
    end

    _show_primer_common(io, primer)
end

function Base.show(io::IO, ::MIME"text/plain", pp::Pair{<:AbstractPrimer})
    fwd, rev = pp.first, pp.second

    # Validate that this is a proper forward/reverse pair
    if !is_forward_primer(fwd) || is_forward_primer(rev)
        invoke(show, Tuple{IO, MIME"text/plain", Pair}, io, MIME"text/plain"(), pp)
        return
    end

    try
        msa = msa_reference(fwd)
        if msa_reference(rev) !== msa
            invoke(show, Tuple{IO, MIME"text/plain", Pair}, io, MIME"text/plain"(), pp)
            return
        end

        N = nseqs(msa)
        L = length(msa)
        amp_start = position_range(fwd).start
        amp_end = position_range(rev).stop
        amp_len = amp_end - amp_start + 1
        overlap = position_range(fwd).stop >= position_range(rev).start ? "!!! OVERLAPPING !!! " : ""
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
        label = position_range(fwd).stop >= position_range(rev).start ? "" : "$(amp_len)bp"
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
        println(io, "Forward: ", consensus_sequence(fwd), " at ", position_range(fwd).start, ":", position_range(fwd).stop)
        println(io, "Reverse: ", consensus_sequence(rev), " at ", position_range(rev).start, ":", position_range(rev).stop)
        mean_tm = (melting_temperature(fwd).mean + melting_temperature(rev).mean) / 2
        delta_tm = abs(melting_temperature(fwd).mean - melting_temperature(rev).mean) / 2
        print(io, "Tm: ", round(mean_tm, digits=1), "±", round(delta_tm, digits=1), " °C")
    catch e
        invoke(show, Tuple{IO, MIME"text/plain", Pair}, io, MIME"text/plain"(), pp)
    end
end
