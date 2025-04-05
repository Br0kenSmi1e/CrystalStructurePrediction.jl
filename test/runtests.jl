using CrystalStructurePrediction
using Test

@testset "interaction.jl" begin
    include("interaction.jl")
end

@testset "build_matrix.jl" begin
    include("build_matrix.jl")
end

@testset "build_problem.jl" begin
    include("build_problem.jl")
end