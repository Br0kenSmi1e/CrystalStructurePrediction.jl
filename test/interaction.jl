using Test
using CrystalStructurePrediction

@testset "lattice" begin
    lattice = Lattice(rand(3,3), [true, true, true])
    @test lattice.vectors ≈ CrystalStructurePrediction.periodic_vectors(lattice)
end

@testset "real_space_sum" begin
    ion_a = Ion("", 1, rand(3))
    ion_b = Ion("", 2, rand(3))
    lattice = Lattice(rand(3,3), [true, true, true])
    alpha = 2.0
    for depth in range(0,4)
        depth_list = [depth for _ in range(1,3)]
        @test real_space_sum(depth_list, ion_a, ion_b, lattice, alpha) ≈ real_space_sum(depth_list, ion_b, ion_a, lattice, alpha)
    end
end

@testset "reciprocal_space_sum" begin
    ion_a = Ion("", 1, rand(3))
    ion_b = Ion("", 2, rand(3))
    lattice = Lattice(rand(3,3), [true, true, true])
    alpha = 2.0
    for depth in range(0,4)
        depth_list = [depth for _ in range(1,3)]
        @test reciprocal_space_sum(depth_list, ion_a, ion_b, lattice, alpha) ≈ reciprocal_space_sum(depth_list, ion_b, ion_a, lattice, alpha)
    end
end

@testset "buckingham_sum" begin
    ion_a = Ion("Sr", 1, rand(3))
    ion_b = Ion("O", 2, rand(3))
    lattice = Lattice(rand(3,3), [true, true, true])
    for depth in range(0,4)
        depth_list = [depth for _ in range(1,3)]
        @test buckingham_sum(depth_list, ion_a, ion_b, lattice) ≈ buckingham_sum(depth_list, ion_b, ion_a, lattice)
    end
end