using Test
using LinearAlgebra: det
using CrystalStructurePrediction, StaticArrays
using CrystalStructurePrediction: interaction_energy

@testset "lattice" begin
    lattice = Lattice(rand(3,3), (true, true, true))
    @test lattice isa Lattice
end

@testset "interaction_energy" begin
    lattice = Lattice([1.0 0.0 0.0; 0.0 0.8 0.0; 0.0 0.6 1.2], (true, true, true))
    ion_a = Ion(IonType(:O, -2, 1.35), [1.0, 0.5, 0.5])
    ion_b = Ion(IonType(:O, -2, 1.35), [0.0, 0.0, 0.2])
    alpha = 2.0 / (abs(det(lattice.vectors)))^(1/3)
    depth = (4, 4, 4)
    @test interaction_energy(ion_a, ion_b, lattice, alpha, depth, depth, depth) ≈ interaction_energy(ion_b, ion_a, lattice, alpha, depth, depth, depth)
    @test interaction_energy(ion_a, ion_b, lattice, alpha, depth, depth, depth) ≈ 322.4479207721009
end

