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
        println(io, "(none)")
    else
        println(io, "\"", description(olig), "\"")
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
        println(io, "(none)")
    else
        println(io, "\"", description(deg), "\"")
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
        println(io, "(none)")
    else
        println(io, "\"", description(go), "\"")
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
    println(io, "  Parent description: ", description(ov))
end