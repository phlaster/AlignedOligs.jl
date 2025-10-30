_truncate_seq(seq::AbstractString, max_width::Int=20) = length(seq) > max_width ?
    seq[1:max(0, max_width-3)] * "..." : seq

function _show_header(io::IO, olig::AbstractOlig)
    println(io, typeof(olig))
    println(io, "  Sequence: ", String(olig))
    println(io, "  Length: ", length(olig))
end

function _show_common_fields(io::IO, olig::AbstractOlig)
    print(io, "  Description: ")
    if isempty(description(olig))
        print(io, "(none)")
    else
        print(io, "\"", description(olig), "\"")
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
    _show_header(io, olig)
    _show_common_fields(io, olig)
end

function Base.show(io::IO, deg::DegenOlig)
    seq_display = _truncate_seq(String(deg))
    print(io, "DegenOlig(\"", seq_display, "\", len=", length(deg))
    print(io, ", n_deg=", n_deg_pos(deg))
    vars_str = n_unique_oligs(deg) > 10000 ? ">10k" : string(n_unique_oligs(deg))
    print(io, ", vars=", vars_str)
    if !isempty(description(deg))
        print(io, ", desc=\"", description(deg), "\"")
    end
    print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", deg::DegenOlig)
    _show_header(io, deg)
    println(io, "  Degenerate positions: ", n_deg_pos(deg))
    println(io, "  Unique variants: ", n_unique_oligs(deg))
    _show_common_fields(io, deg)
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
    println(io, "  Gaps: ", length(go.gaps))
    _show_common_fields(io, go)
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
    println(io, typeof(ov))
    println(io, "  Viewed sequence: ", String(ov))
    println(io, "  Length: ", length(ov))
    println(io, "  Range: ", ov.range)
    println(io, "  Parent description: ", description(ov))
end