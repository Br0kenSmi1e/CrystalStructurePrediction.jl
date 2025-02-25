using CrystalStructurePrediction
using Test

# @testset "CrystalStructurePrediction.jl" begin
    # Write your tests here.
# end

@testset "interaction.jl" begin
    include("interaction.jl")
end

@testset "build_matrix.jl" begin
    include("build_matrix.jl")
end
