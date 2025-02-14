using CrystalStructurePrediction
using Documenter

DocMeta.setdocmeta!(CrystalStructurePrediction, :DocTestSetup, :(using CrystalStructurePrediction); recursive=true)

makedocs(;
    modules=[CrystalStructurePrediction],
    authors="Br0kenSmi1e",
    sitename="CrystalStructurePrediction.jl",
    format=Documenter.HTML(;
        canonical="https://Br0kenSmi1e.github.io/CrystalStructurePrediction.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Br0kenSmi1e/CrystalStructurePrediction.jl",
    devbranch="main",
)
