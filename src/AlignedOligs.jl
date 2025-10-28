module AlignedOligs

"""
    Package AlignedOligs

Nucleic acid oligomers aligning and PCR primers construction

$(isnothing(get(ENV, "CI", nothing)) ? ("\nPackage local path: $(pathof(AlignedOligs))") : "") 
"""
AlignedOligs

using Reexport
using FastaIO
using Statistics
using Random
using Statistics
using ProgressMeter

@reexport using SeqFold
using SeqFold: Structure

export AbstractOlig, Olig, DegenerateOlig, OligView, GappedOlig
export n_unique_oligs, n_deg_pos, nondegens, description, hasgaps
export AbstractMSA, MSA, MSAView, msadepth, msadet, width, height, root, bval#, min_tresh
export nseqs, width, getsequence, getsequence, get_base_count, consensus_major, consensus_degen, dry_msa
export Primer, construct_primers, best_pairs, nucleotide_diversity

include("utils.jl")
include("Oligs.jl")
include("MSA.jl")
include("SeqFold_methods.jl")
include("Primers.jl")
include("show.jl")

end
