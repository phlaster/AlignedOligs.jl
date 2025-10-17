const IUPAC_B2V = Dict(
    'A'=>Tuple("A"),  'C'=>Tuple("C"),  'G'=>Tuple("G"),  'T'=>Tuple("T"),
    'R'=>Tuple("AG"), 'Y'=>Tuple("CT"), 'S'=>Tuple("CG"), 'W'=>Tuple("AT"),
    'K'=>Tuple("GT"), 'M'=>Tuple("AC"), 'B'=>Tuple("CGT"),'D'=>Tuple("AGT"),
    'H'=>Tuple("ACT"),'V'=>Tuple("ACG"),'N'=>Tuple("ACGT")
)
const IUPAC_V2B = Dict(
    collect("A")=>'A', collect("C")=>'C', collect("G")=>'G', collect("T")=>'T',
    collect("AC")=>'M', collect("AG")=>'R', collect("AT")=>'W', collect("CG")=>'S',
    collect("CT")=>'Y', collect("GT")=>'K', collect("ACG")=>'V', collect("ACT")=>'H',
    collect("AGT")=>'D', collect("CGT")=>'B', collect("ACGT")=>'N',
)

const IUPAC_COUNTS = Dict(
    'A'=>1, 'C'=>1, 'G'=>1, 'T'=>1,
    'M'=>2, 'R'=>2, 'W'=>2, 'S'=>2, 'Y'=>2, 'K'=>2,
    'V'=>3, 'H'=>3, 'D'=>3, 'B'=>3, 'N'=>4
)

const IUPAC_GC_CONTENT = Dict(
    'A'=>0.0, 'T'=>0.0, 'C'=>1.0, 'G'=>1.0,
    'R'=>0.5, 'Y'=>0.5, 'S'=>1.0, 'W'=>0.0,
    'K'=>0.5, 'M'=>0.5, 'B'=>2/3, 'D'=>1/3,
    'H'=>1/3, 'V'=>2/3, 'N'=>0.5
)

const NON_DEGEN_BASES = Tuple("ACGT")
const DEGEN_BASES = Tuple("MRWSYKVHDBN")
const ALL_BASES = Tuple("ACGTMRWSYKVHDBN")

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
    d = Dict{Char, Vector{Float64}}()
    bases = "ACGT"
    for (k, v) in IUPAC_B2V
        probs = zeros(Float64, 4)
        for base in v
            idx = findfirst(==(base), bases)
            probs[idx] += 1.0
        end
        if !isempty(v)
            probs ./= length(v)
        end
        d[k] = probs
    end
    d
end