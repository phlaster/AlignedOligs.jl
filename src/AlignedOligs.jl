module AlignedOligs

"""
    Package AlignedOligs

Nucleic acid oligomers aligning and PCR primers construction

$(isnothing(get(ENV, "CI", nothing)) ? ("\nPackage local path: $(pathof(AlignedOligs))") : "") 
"""
AlignedOligs

using Reexport
using FastaIO
using MAFFT_jll
using Statistics
using Random

@reexport using SeqFold
using SeqFold: Structure

export AbstractOlig, Olig, DegenerateOlig, OligView, n_unique_oligs, n_deg_pos, nondegens, description
# export MSA, nseqs, getrange, getsequence, get_base_count, dry_msa, consensus_major, consensus_degen

include("utils.jl")
include("Oligs.jl")
# include("MSA.jl")



end
