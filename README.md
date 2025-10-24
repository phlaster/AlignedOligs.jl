# AlignedOligs

[![Build Status](https://github.com/phlaster/AlignedOligs.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/phlaster/AlignedOligs.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

# Quickstart

```julia
using AlignedOligs, MAFFT_jll

nal = "ActR.fasta"

m = MSA(nal, mafft=true, minor_threshold=0.05, bootstrap=200);

pfs = construct_primers(m);
prs = construct_primers(m, is_forward=false);
amplicons = best_pairs(pfs, prs, amplicon_len=150:999);

open("out.txt", "w") do f
    for pp in amplicons
        show(f, MIME"text/plain"(), pp)
        println(f)
        println(f)
    end
end
```