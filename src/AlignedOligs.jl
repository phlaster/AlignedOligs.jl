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
using Statistics
using ProgressMeter

@reexport using SeqFold
using SeqFold: Structure

export AbstractOlig, Olig, DegenerateOlig, OligView, GappedOlig
export n_unique_oligs, n_deg_pos, nondegens, description, hasgaps
export AbstractMSA, MSA, MSAView
export nseqs, width, getsequence, getsequence, get_base_count, consensus_major, consensus_degen, dry_msa

include("utils.jl")
include("Oligs.jl")
include("Oligs_display.jl")
include("MSA.jl")
include("MSA_display.jl")
include("SeqFold_methods.jl")



end
