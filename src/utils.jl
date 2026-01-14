const NON_DEGEN_BASES = "ACGT" |> collect |> sort |> Tuple
const DEGEN_BASES = "MRWSYKVHDBN" |> collect |> sort |> Tuple
const ALL_BASES = [NON_DEGEN_BASES..., DEGEN_BASES...] |> sort |> Tuple
const BASES_W_GAPS = (ALL_BASES..., '-')

const IUPAC_B2V = Dict(
    'A'=>Tuple("A"),  'C'=>Tuple("C"),  'G'=>Tuple("G"),  'T'=>Tuple("T"),
    'R'=>Tuple("AG"), 'Y'=>Tuple("CT"), 'S'=>Tuple("CG"), 'W'=>Tuple("AT"),
    'K'=>Tuple("GT"), 'M'=>Tuple("AC"), 'B'=>Tuple("CGT"),'D'=>Tuple("AGT"),
    'H'=>Tuple("ACT"),'V'=>Tuple("ACG"),'N'=>Tuple("ACGT")
)
const IUPAC_V2B = Dict(v=>k for (k,v) in IUPAC_B2V)
const IUPAC_COUNTS = Dict(k=>length(v) for (k,v) in IUPAC_B2V)
const IUPAC_GC_CONTENT = Dict(k=>count(in("CG"), v)/length(v) for (k,v) in IUPAC_B2V)

const DNA_COMP_TABLE_DEG = let 
    _comp_dna_deg(b::UInt8)::UInt8 =
    b == UInt8('A') ? UInt8('T') : b == UInt8('a') ? UInt8('T') :
    b == UInt8('T') ? UInt8('A') : b == UInt8('t') ? UInt8('A') :
    b == UInt8('C') ? UInt8('G') : b == UInt8('c') ? UInt8('G') :
    b == UInt8('G') ? UInt8('C') : b == UInt8('g') ? UInt8('C') :
    b == UInt8('R') ? UInt8('Y') : b == UInt8('r') ? UInt8('Y') :
    b == UInt8('Y') ? UInt8('R') : b == UInt8('y') ? UInt8('R') :
    b == UInt8('S') ? UInt8('S') : b == UInt8('s') ? UInt8('S') :
    b == UInt8('W') ? UInt8('W') : b == UInt8('w') ? UInt8('W') :
    b == UInt8('K') ? UInt8('M') : b == UInt8('k') ? UInt8('M') :
    b == UInt8('M') ? UInt8('K') : b == UInt8('m') ? UInt8('K') :
    b == UInt8('B') ? UInt8('V') : b == UInt8('b') ? UInt8('V') :
    b == UInt8('D') ? UInt8('H') : b == UInt8('d') ? UInt8('H') :
    b == UInt8('H') ? UInt8('D') : b == UInt8('h') ? UInt8('D') :
    b == UInt8('V') ? UInt8('B') : b == UInt8('v') ? UInt8('B') :
    b == UInt8('N') ? UInt8('N') : b == UInt8('n') ? UInt8('N') :
    UInt8('-')
    [_comp_dna_deg(UInt8(i)) for i in 0:255]
end

const IUPAC_PROBS = let
    d = Dict{Char, NTuple{4, Float64}}()
    for (k, v) in IUPAC_B2V
        probs = zeros(Float64, 4)
        for base in v
            idx = findfirst(==(base), NON_DEGEN_BASES)
            probs[idx] += 1.0
        end
        if !isempty(v)
            probs ./= length(v)
        end
        d[k] = Tuple(probs)
    end
    d
end

const MAX_GC_OPTIONS = let
    d = Dict{Char, Tuple}()
    for (k0, v0) in IUPAC_B2V
        _v = filter(in("GC"), v0)
        v = isempty(_v) ? v0 : _v
        d[k0] = v
    end
    d
end

const MIN_GC_OPTIONS = let
    d = Dict{Char, Tuple}()
    for (k0, v0) in IUPAC_B2V
        _v = filter(!in("GC"), v0)
        v = isempty(_v) ? v0 : _v
        d[k0] = v
    end
    d
end