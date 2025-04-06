using CrystalStructurePrediction
using Test, Documenter

@testset "struct.jl" begin
    include("struct.jl")
end

@testset "interaction.jl" begin
    include("interaction.jl")
end

@testset "build_problem.jl" begin
    include("build_problem.jl")
end

# Run doctests
DocMeta.setdocmeta!(CrystalStructurePrediction, :DocTestSetup, :(using CrystalStructurePrediction); recursive=true)
doctest(CrystalStructurePrediction)