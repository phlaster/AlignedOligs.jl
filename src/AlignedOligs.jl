"""
    AlignedOligs.jl

Package for nucleic acid oligomers aligning and PCR primers construction.


$(isnothing(get(ENV, "CI", nothing)) ? ("Package local path: " * string(pathof(AlignedOligs))) : "")
"""
module AlignedOligs

using Reexport
# using FastaIO
# using Statistics
# using Statistics
# using ProgressMeter

@reexport using SeqFold

export Oligs, Primers, Alignments

# export AbstractPrimer

# export Primer
# export construct_primers, best_pairs

include("Oligs.jl")
include("MSA.jl")
# include("Primers.jl")
# include("SeqFold_methods.jl")
# include("show.jl")

end # module