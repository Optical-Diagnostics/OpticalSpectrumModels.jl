using SpectraUtils
using Documenter

DocMeta.setdocmeta!(SpectraUtils, :DocTestSetup, :(using SpectraUtils); recursive=true)

makedocs(;
    modules=[SpectraUtils],
    authors="Jan Kuhfeld <jankuhfeld@plasma-matters.nl> and contributors",
    sitename="SpectraUtils.jl",
    format=Documenter.HTML(;
        canonical="https://optical-diagnostics.github.io/SpectraUtils.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Spectra" => "spectra.md",
        "Linewidth Estimates" => "linewidths.md",
        "HITRAN" => "hitran.md",
    ],
)

deploydocs(;
    repo="github.com/Optical-Diagnostics/SpectraUtils.jl",
    devbranch="main",
)
