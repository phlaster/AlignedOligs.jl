using AlignedOligs
using Documenter

DocMeta.setdocmeta!(AlignedOligs, :DocTestSetup, :(using AlignedOligs); recursive=true)

makedocs(;
    modules=[AlignedOligs],
    authors="phlaster <phlaster@users.noreply.github.com>",
    sitename="AlignedOligs.jl",
    format=Documenter.HTML(;
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
