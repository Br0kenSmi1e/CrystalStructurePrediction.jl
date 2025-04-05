using Test
using CrystalStructurePrediction, StaticArrays
using CrystalStructurePrediction: interaction_energy, real_space_Ewald_sum, reciprocal_space_Ewald_sum, buckingham_sum

@testset "lattice" begin
    lattice = Lattice(rand(3,3), (true, true, true))
    @test lattice isa Lattice
end

@testset "real_space_sum" begin
    ion_a = Ion(IonType(:none, 1, 0.0), [0.0, 0.25, 0.0])
    ion_b = Ion(IonType(:none, 2, 0.0), [0.5, 0.0, 0.25])
    lattice = Lattice([1.1 0.0 0.0; 0.0 0.5 0.9; 0.0 1.0 0.8], (true, true, true))
    alpha = 2.0
    expected = [-16.24825982870472, -9.546793813049247, -9.353356845659718, -9.351227627445452, -9.351217663090333]
    for depth in range(0,4)
        depth_list = ntuple(x->depth, 3)
        @test real_space_Ewald_sum(ion_a, ion_b, lattice, alpha, depth_list) ≈ real_space_Ewald_sum(ion_b, ion_a, lattice, alpha, depth_list)
        @test real_space_Ewald_sum(ion_a, ion_a, lattice, alpha, depth_list) ≈ expected[depth + 1]
    end
end

@testset "reciprocal_space_sum" begin
    ion_a = Ion(IonType(:none, 1, 0.0), [0.0, 0.25, 0.0])
    ion_b = Ion(IonType(:none, 2, 0.0), [0.5, 0.0, 0.25])
    lattice = Lattice([1.1 0.0 0.0; 0.0 0.5 0.9; 0.0 1.0 0.8], (true, true, true))
    alpha = 2.0
    expected = [0.0, 2.944589381733499, 2.946645781289601, 2.9466458258563697, 2.946645825856462]
    for depth in range(0,4)
        depth_list = ntuple(x->depth, 3)
        @test reciprocal_space_Ewald_sum(ion_a, ion_b, lattice, alpha, depth_list) ≈ reciprocal_space_Ewald_sum(ion_b, ion_a, lattice, alpha, depth_list)
        @test reciprocal_space_Ewald_sum(ion_a, ion_a, lattice, alpha, depth_list) ≈ expected[depth + 1]
    end
end

@testset "buckingham_sum" begin
    lattice = Lattice([1.0 0.0 0.0; 0.0 0.8 0.0; 0.0 0.6 1.2], (true, true, true))
    ion_a = Ion(IonType(:Sr, 1, 0.0), [0.0, 0.25, 0.0])
    ion_b = Ion(IonType(:O, 2, 0.0), [0.5, 0.0, 0.25])
    lattice = Lattice([1.1 0.0 0.0; 0.0 0.5 0.9; 0.0 1.0 0.8], (true, true, true))
    expected = zeros(5)
    for depth in range(0,4)
        depth_list = ntuple(x->depth, 3)
        @test buckingham_sum(ion_a, ion_b, lattice, depth_list) ≈ buckingham_sum(ion_b, ion_a, lattice, depth_list)
        @test buckingham_sum(ion_a, ion_a, lattice, depth_list) ≈ expected[depth + 1]
    end
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

