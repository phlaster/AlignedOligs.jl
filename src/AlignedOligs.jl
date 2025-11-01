module AlignedOligs

"""
    Package AlignedOligs

Nucleic acid oligomers aligning and PCR primers construction

$(isnothing(get(ENV, "CI", nothing)) ? ("\nPackage local path: $(pathof(AlignedOligs))") : "") 
"""
AlignedOligs

using Reexport

@reexport using SeqFold

export Oligs, Primers, Alignments


include("Oligs.jl")
include("MSA.jl")
include("Primers.jl")

end # module