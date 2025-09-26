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
    if !isempty(olig.description)
        print(io, ", desc=\"", olig.description, "\"")
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
    print(io, ", n_deg=", deg.n_deg_pos)
    vars_str = deg.n_unique_oligs > 10000 ? ">10k" : string(deg.n_unique_oligs)
    print(io, ", vars=", vars_str)
    if !isempty(deg.description)
        print(io, ", desc=\"", deg.description, "\"")
    end
    print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", deg::DegenerateOlig)
    println(io, "DegenerateOlig")
    println(io, "  Sequence: ", String(deg))
    println(io, "  Length: ", length(deg))
    println(io, "  Degenerate positions: ", deg.n_deg_pos)
    println(io, "  Unique variants: ", deg.n_unique_oligs)
    print(io, "  Description: ")
    if isempty(deg.description)
        println(io, "(none)")
    else
        println(io, "\"", deg.description, "\"")
    end
end

function Base.show(io::IO, go::GappedOlig)
    seq_display = _truncate_seq(String(go))
    print(io, "GappedOlig(\"", seq_display, "\", len=", length(go))
    print(io, ", gaps=", length(go.gaps))
    if !isempty(go.description)
        print(io, ", desc=\"", go.description, "\"")
    end
    print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", go::GappedOlig)
    println(io, "GappedOlig")
    println(io, "  Gapped sequence: ", String(go))
    println(io, "  Length (with gaps): ", length(go))
    println(io, "  Underlying Olig: ", go.olig)
    println(io, "  Gaps: ", go.gaps)
    print(io, "  Description: ")
    if isempty(go.description)
        println(io, "(none)")
    else
        println(io, "\"", go.description, "\"")
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