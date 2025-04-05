using CrystalStructurePrediction, Test

@testset "build_problem" begin
    # Crystal structure parameters
    grid_size = (2, 2, 2)
    population_list = [1, 1, 3]  # 1 Sr, 1 Ti, 3 O atoms
    Sr, Ti, O = IonType(:Sr, 2, 1.18), IonType(:Ti, 4, 0.42), IonType(:O, -2, 1.35)
    type_list = [Sr, Ti, O]
    
    # Lattice parameters
    lattice_constant = 3.899  # Å
    L = lattice_constant * [1 0 0; 0 1 0; 0 0 1]
    lattice = Lattice(L, (true, true, true))
    
    # Ewald summation parameters
    depth = (4, 4, 4)
    alpha = 2 / lattice_constant
    
    # Build ion list and proximal pairs
    ion_list = ions_on_grid(grid_size, type_list)
    proximal_pairs = build_proximal_pairs(ion_list, lattice, 0.75)
    
    # Solve the linear problem
    res_linear = build_linear_problem(ion_list, Dict(type_list .=> population_list), lattice) do ion_a, ion_b, lattice
        interaction_energy(ion_a, ion_b, lattice, alpha, depth, depth, depth)
    end
    @test res_linear.energy ≈ -6.061349350569214

    # The quadratic problem formulation
    res_quadratic = build_quadratic_problem(ion_list, Dict(type_list .=> population_list), lattice) do ion_a, ion_b, lattice
        interaction_energy(ion_a, ion_b, lattice, alpha, depth, depth, depth)
    end
    @test res_quadratic.energy ≈ -6.061349350569214/2
end