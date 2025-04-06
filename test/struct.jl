using CrystalStructurePrediction, Test
using CrystalStructurePrediction: species, charge, radii, fractional, cartesian, minimum_distance, distance

@testset "struct.jl" begin
    # Test IonType constructor
    ion_type = IonType(:O, -2, 1.35)
    @test species(ion_type) == :O
    @test charge(ion_type) == -2
    @test radii(ion_type) == 1.35
    
    # Test Ion constructor
    ion = Ion(ion_type, [1//2, 1//2, 1//2])
    @test species(ion) == :O
    @test charge(ion) == -2
    @test radii(ion) == 1.35
    @test ion.frac_pos ≈ [1//2, 1//2, 1//2]
    
    # Test Lattice constructor
    vectors = [0.5 0 0; 0 0.5 0; 0 0 0.5]
    lattice = Lattice(vectors)
    @test lattice.vectors ≈ vectors
    @test minimum_distance([0//1, 0, 0], [9//10, 9//10, 9//10], lattice) ≈ 0.05 * sqrt(3)
    @test distance([0//1, 0, 0], [9//10, 9//10, 9//10]) ≈ sqrt(3) * 0.9
    @test distance(lattice, [0//1, 0, 0], [9//10, 9//10, 9//10]) ≈ sqrt(3) * 0.45

    # Test fractional and cartesian conversion
    vectors = randn(3, 3)
    lattice = Lattice(vectors)
    frac_pos = rand(3)
    cart_pos = cartesian(lattice, frac_pos)
    @test fractional(lattice, cart_pos) ≈ frac_pos
end