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
    max_display_height = max(5, terminal_height - 6)
    max_display_width = max(20, terminal_width - 30)

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
                println(io, "...")
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
                    println(io, start_part * "..." * end_part)
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
                     println(io, start_part * "..." * end_part)
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
                    println(io, start_part * "..." * end_part)
                end
            end
        end
    end
end