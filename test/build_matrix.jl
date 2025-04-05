using Test
using LinearAlgebra: det
using CrystalStructurePrediction
using CrystalStructurePrediction: interaction_energy, build_matrix, build_ion_list

@testset "interaction_energy" begin
    lattice = Lattice(rand(3,3), (true, true, true))
    ion_a = Ion(IonType(:O, -2, 1.35), rand(3))
    ion_b = Ion(IonType(:O, -2, 1.35), rand(3))
    alpha = 2.0 / (abs(det(lattice.vectors)))^(1/3)
    depth = (4, 4, 4)
    @test interaction_energy(ion_a, ion_b, lattice, alpha, depth, depth, depth) ≈ interaction_energy(ion_b, ion_a, lattice, alpha, depth, depth, depth)
end

@testset "build_matrix" begin
    lattice = Lattice(rand(2, 2), (true, true))
    ion_list = build_ion_list((1, 1), [IonType(:None, 1, 0.0)])
    alpha = 2.0 / (abs(det(lattice.vectors)))^(1/2)
    depth = (0, 0)
    @test build_matrix(ion_list, lattice, interaction_energy, (alpha, depth, depth, depth)) ≈ [-14.399645351950543*alpha/(π)^0.5]
end