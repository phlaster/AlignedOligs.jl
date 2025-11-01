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

const histbars = ['▁', '▂', '▃', '▄', '▅', '▆', '▇', '█']

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
    max_display_height = max(5, terminal_height - 9)  # Increased by 2 to account for hist bars and number line
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

    # Add histogram bars for depth with colored major nucleotide
    if has_desc
        print(io, " "^desc_width)  # Print spaces to align with descriptions
    end
    
    hist_line = IOBuffer()
    for dj in 1:length(displayed_cols_range)
        # Get major nucleotide and depth for this column
        pos_counts = zeros(Float64, 4)
        for i in 1:n_sequences
            c = getsequence(msa, i, dj)
            probs = IUPAC_PROBS[c]
            pos_counts .+= probs
        end
        
        # Find the major nucleotide (highest probability base)
        max_prob, max_idx = findmax(pos_counts)
        major_nuc = ['A', 'C', 'G', 'T'][max_idx]
        depth = sum(pos_counts)
        
        # Normalize depth to 0-1 scale for bar selection (assuming max possible depth is n_sequences)
        max_possible_depth = n_sequences
        normalized_depth = depth / max_possible_depth
        
        # Select bar character based on normalized depth (0-1) -> index 1-8
        bar_index = clamp(floor(Int, normalized_depth * 8) + 1, 1, 8)
        bar_char = histbars[bar_index]
        
        # Get color for the major nucleotide
        color = get(BASE_COLORS, major_nuc, :normal)
        
        # Print the bar character with the color of the major nucleotide
        printstyled(io, bar_char; color=color)
    end
    println(io)

    # Print sequence lines
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

    # Add column number line
    if has_desc
        print(io, " "^desc_width)  # Print spaces to align with descriptions
    end
    
    # Calculate start and end absolute positions
    start_abs = abs_cols[1]
    end_abs = abs_cols[end]
    
    # Create a character array for the number line, same length as displayed sequence
    num_line_chars = fill(' ', length(displayed_cols_range))
    
    # Place the start number - first digit at first position
    start_num_str = string(start_abs)
    start_num_len = length(start_num_str)
    if length(num_line_chars) >= start_num_len
        for i in 1:start_num_len
            if i <= length(num_line_chars)
                num_line_chars[i] = start_num_str[i]
            end
        end
    end
    
    # Place the end number - last digit at last position
    end_num_str = string(end_abs)
    end_num_len = length(end_num_str)
    if length(num_line_chars) >= end_num_len
        for i in 1:end_num_len
            idx = length(num_line_chars) - end_num_len + i
            if idx >= 1
                num_line_chars[idx] = end_num_str[i]
            end
        end
    end
    
    # Place '*' symbols at positions where the absolute position is divisible by 10
    # and the character at that position is still a space (not overwritten by numbers)
    for (rel_idx, abs_pos) in enumerate(abs_cols)
        if abs_pos % 10 == 0 && num_line_chars[rel_idx] == ' '
            num_line_chars[rel_idx] = '*'
        end
    end
    
    print(io, String(num_line_chars))
    println(io)

    if n_display_seqs < n_sequences
        println(io, "...")
    end
end

function _is_full_height(msa::AbstractMSA)
    depth = msadepth(msa, 1)
    return all(d -> d == depth, msadepth(msa, 1:length(msa)))
end