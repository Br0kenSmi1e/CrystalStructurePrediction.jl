using Test
using CrystalStructurePrediction

@testset "real_space_sum" begin
    pos1 = rand(3)
    pos2 = rand(3)
    L = rand(Float64, (3,3))
    alpha = 2.0
    for depth in range(0,8)
        @test real_space_sum(depth, pos1, pos1, L, alpha) ≈ real_space_sum(depth, pos2, pos2, L, alpha)
    end
end