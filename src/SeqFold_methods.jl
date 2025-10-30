
SeqFold.gc_content(primer::AbstractPrimer)::Float64 = SeqFold.gc_content(consensus_oligomer(primer))


SeqFold.dg(primer::AbstractPrimer; kwargs...) = SeqFold.dg(consensus_oligomer(primer); kwargs...)



SeqFold.dg_cache(primer::AbstractPrimer; kwargs...) = SeqFold.dg_cache(consensus_oligomer(primer); kwargs...)



SeqFold.tm(primer1::AbstractPrimer, primer2::AbstractOlig; kwargs...) = 
    SeqFold.tm(consensus_oligomer(primer1), primer2; kwargs...)

SeqFold.tm(olig::AbstractOlig, primer2::AbstractPrimer; kwargs...) = 
    SeqFold.tm(olig, consensus_oligomer(primer2); kwargs...)

SeqFold.tm(primer1::AbstractPrimer, primer2::AbstractPrimer; kwargs...) = 
    SeqFold.tm(consensus_oligomer(primer1), consensus_oligomer(primer2); kwargs...)



SeqFold.tm(primer::AbstractPrimer; kwargs...) = 
    SeqFold.tm(consensus_oligomer(primer); kwargs...)


SeqFold.tm_cache(primer1::AbstractPrimer, primer2::AbstractOlig; kwargs...) = 
    SeqFold.tm_cache(consensus_oligomer(primer1), primer2; kwargs...)

SeqFold.tm_cache(olig::AbstractOlig, primer2::AbstractPrimer; kwargs...) = 
    SeqFold.tm_cache(olig, consensus_oligomer(primer2); kwargs...)

SeqFold.tm_cache(primer1::AbstractPrimer, primer2::AbstractPrimer; kwargs...) = 
    SeqFold.tm_cache(consensus_oligomer(primer1), consensus_oligomer(primer2); kwargs...)



SeqFold.tm_cache(primer::AbstractPrimer; kwargs...) = 
    SeqFold.tm_cache(consensus_oligomer(primer); kwargs...)




SeqFold.dot_bracket(primer::AbstractPrimer, structs::Vector{SeqFold.Structure}) = 
    SeqFold.dot_bracket(consensus_oligomer(primer), structs)



SeqFold.gc_cache(primer::AbstractPrimer)::Matrix{Float64} = 
    SeqFold.gc_cache(consensus_oligomer(primer))


# # ============================================================================
# # GENERIC CONVENIENCE METHODS
# # ============================================================================

# """
#     Base.String(olig::AbstractPrimer)

# Convert primer to string representation of its consensus sequence.
# """
# Base.String(primer::AbstractPrimer) = String(consensus_oligomer(primer))


# """
#     Base.length(primer::AbstractPrimer)

# Get the length of the primer's consensus sequence.
# """
# Base.length(primer::AbstractPrimer) = length(consensus_oligomer(primer))


# """
#     Base.isempty(primer::AbstractPrimer)

# Check if the primer has an empty consensus sequence.
# """
# Base.isempty(primer::AbstractPrimer) = isempty(consensus_oligomer(primer))


# """
#     Base.iterate(primer::AbstractPrimer, state...)

# Iterate over characters in the primer's consensus sequence.
# """
# Base.iterate(primer::AbstractPrimer, state...) = iterate(consensus_oligomer(primer), state...)


# """
#     Base.getindex(primer::AbstractPrimer, i::Int)

# Get character at position i in the primer's consensus sequence.
# """
# Base.getindex(primer::AbstractPrimer, i::Int) = getindex(consensus_oligomer(primer), i)


# """
#     Base.getindex(primer::AbstractPrimer, r::UnitRange{Int})

# Get substring range from the primer's consensus sequence.
# """
# Base.getindex(primer::AbstractPrimer, r::UnitRange{Int}) = getindex(consensus_oligomer(primer), r)


# # ============================================================================
# # TYPE CONVERSION AND COMPATIBILITY
# # ============================================================================

# """
#     Olig(primer::AbstractPrimer)

# Convert a primer to its underlying Olig representation (if non-degenerate).

# Throws: InexactError if primer is degenerate
# """
# function Olig(primer::AbstractPrimer)
#     if n_deg_pos(primer) > 0
#         throw(InexactError(:Olig, Olig, primer))
#     end
#     Olig(String(primer), description(primer))
# end


# """
#     DegenOlig(primer::AbstractPrimer)

# Convert a primer to its DegenOlig representation.
# """
# function DegenOlig(primer::AbstractPrimer)
#     DegenOlig(String(primer), n_deg_pos(primer), n_unique_oligs(primer), description(primer))
# end


# """
#     GappedOlig(primer::AbstractPrimer)

# Convert a primer to a GappedOlig representation.

# Note: Primers don't have gaps by construction, but this provides compatibility.
# """
# GappedOlig(primer::AbstractPrimer) = GappedOlig(DegenOlig(primer), Pair{Int}[])


# """
#     OligView(primer::AbstractPrimer, range::UnitRange{Int})

# Create an OligView of a primer's consensus sequence.
# """
# OligView(primer::AbstractPrimer, range::UnitRange{Int}) = OligView(DegenOlig(primer), range)


# # ============================================================================
# # SPECIALIZED METHODS FOR EXTENSIBILITY
# # ============================================================================

# """
#     has_thermodynamic_data(olig)

# Check if an oligonucleotide has pre-calculated thermodynamic properties.

# Supports: AbstractOlig (always false), AbstractPrimer (always true via stored data)
# """
# has_thermodynamic_data(::AbstractOlig) = false
# has_thermodynamic_data(::AbstractPrimer) = true


# """
#     get_tm_stats(olig)

# Get melting temperature statistics if available.

# Supports: AbstractPrimer (returns stored data), AbstractOlig (calculates on-demand)
# """
# get_tm_stats(olig::AbstractOlig; kwargs...) = SeqFold.tm(olig; kwargs...)
# get_tm_stats(primer::AbstractPrimer) = melting_temperature(primer)


# """
#     get_dg_value(olig)

# Get free energy value if available.

# Supports: AbstractPrimer (returns stored data), AbstractOlig (calculates on-demand)
# """
# get_dg_value(olig::AbstractOlig; kwargs...) = SeqFold.dg(olig; kwargs...)
# get_dg_value(primer::AbstractPrimer) = free_energy(primer)


# """
#     get_gc_stats(olig)

# Get GC content value if available.

# Supports: AbstractPrimer (returns stored data), AbstractOlig (calculates on-demand)
# """
# get_gc_stats(olig::AbstractOlig) = SeqFold.gc_content(olig)
# get_gc_stats(primer::AbstractPrimer) = gc_content(primer)