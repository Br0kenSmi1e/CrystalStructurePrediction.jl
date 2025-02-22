using Test
using CrystalStructurePrediction

@testset "lattice" begin
    lattice = Lattice(rand(3,3), [true, true, true])
    @test lattice.vectors ≈ CrystalStructurePrediction.periodic_vectors(lattice)
end

@testset "real_space_sum" begin
    pos1 = rand(3)
    pos2 = rand(3)
    lattice = Lattice(rand(3,3), [true, true, true])
    alpha = 2.0
    for depth in range(0,4)
        depth_list = [depth for _ in range(1,3)]
        @test real_space_sum(depth_list, pos1, pos1, lattice, alpha) ≈ real_space_sum(depth_list, pos2, pos2, lattice, alpha)
    end
end