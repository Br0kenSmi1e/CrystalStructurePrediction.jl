using CrystalStructurePrediction, Test
using CrystalStructurePrediction: species, charge, radii

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
    vectors = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    pbc = (true, true, true)
    lattice = Lattice(vectors, pbc)
    @test lattice.vectors ≈ vectors
    @test lattice.pbc == pbc
end