using CrystalStructurePrediction, Test
using CrystalStructurePrediction: species, charge, radii, fractional, cartesian, minimum_distance, distance

@testset "struct.jl" begin
    # Test IonType constructor
    ion_type = IonType(:O, -2, 1.35)
    @test species(ion_type) == :O
    @test charge(ion_type) == -2
    @test radii(ion_type) == 1.35
    
    # Test Ion constructor
    ion = Ion(ion_type, [0.5, 0.5, 0.5])
    @test species(ion) == :O
    @test charge(ion) == -2
    @test radii(ion) == 1.35
    @test ion.frac_pos ≈ [0.5, 0.5, 0.5]
    
    # Test Lattice constructor
    vectors = [0.5 0.0 0.0; 0.0 0.5 0.0; 0.0 0.0 0.5]
    pbc = (true, true, true)
    lattice = Lattice(vectors, pbc)
    @test lattice.vectors ≈ vectors
    @test lattice.pbc == pbc
    @test minimum_distance([0.0, 0.0, 0.0], [0.9, 0.9, 0.9], lattice) ≈ 0.05 * sqrt(3)
    @test distance([0.0, 0.0, 0.0], [0.9, 0.9, 0.9]) ≈ sqrt(3) * 0.9
    @test distance(lattice, [0.0, 0.0, 0.0], [0.9, 0.9, 0.9]) ≈ sqrt(3) * 0.45

    # Test fractional and cartesian conversion
    vectors = randn(3, 3)
    lattice = Lattice(vectors, pbc)
    frac_pos = rand(3)
    cart_pos = cartesian(lattice, frac_pos)
    @test fractional(lattice, cart_pos) ≈ frac_pos
end