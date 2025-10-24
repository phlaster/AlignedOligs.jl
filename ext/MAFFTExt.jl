__precompile__(false)
module MAFFTExt

using AlignedOligs
using MAFFT_jll
using FastaIO

function AlignedOligs._align!(fasta_content::Vector{Tuple{String,String}})
    buffer = IOBuffer()
    writefasta(buffer, fasta_content; check_description=false)
    frombuffer = take!(buffer)

    proc = open(`$(MAFFT_jll.mafft()) --quiet -`, "r+")
    # proc = open(pipeline(`$(MAFFT_jll.mafft()) -`, stderr=devnull), "r+")
    write(proc, frombuffer)
    close(proc.in)
    aligned_output = read(proc, String)

    fasta_content .= readfasta(IOBuffer(aligned_output))
    return
end

end  # module
