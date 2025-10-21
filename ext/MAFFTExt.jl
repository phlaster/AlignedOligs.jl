__precompile__(false)
module MAFFTExt

using AlignedOligs
using MAFFT_jll
using FastaIO

function AlignedOligs._align_if_needed!(seqs_desc::Vector{Tuple{String,String}}, mafft::Bool)
    mafft || return

    buffer = IOBuffer()
    writefasta(buffer, seqs_desc; check_description=false)
    fasta_content = take!(buffer)

    proc = open(`$(MAFFT_jll.mafft()) --quiet -`, "r+")
    # proc = open(pipeline(`$(MAFFT_jll.mafft()) -`, stderr=devnull), "r+")
    write(proc, fasta_content)
    close(proc.in)
    aligned_output = read(proc, String)

    seqs_desc .= readfasta(IOBuffer(aligned_output))
    return
end

end  # module MAFFTExt