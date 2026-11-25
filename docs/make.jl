using OpticalSpectrumModels
using Documenter

DocMeta.setdocmeta!(OpticalSpectrumModels, :DocTestSetup, :(using OpticalSpectrumModels); recursive=true)

makedocs(;
    modules=[OpticalSpectrumModels],
    authors="Jan Kuhfeld <jankuhfeld@plasma-matters.nl> and contributors",
    sitename="OpticalSpectrumModels.jl",
    format=Documenter.HTML(;
        canonical="https://optical-diagnostics.github.io/OpticalSpectrumModels.jl",
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
    repo="github.com/Optical-Diagnostics/OpticalSpectrumModels.jl",
    devbranch="main",
)
