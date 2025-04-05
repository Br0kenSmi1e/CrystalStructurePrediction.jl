using CrystalStructurePrediction, Test

@testset "build_problem" begin
    # Crystal structure parameters
    grid_size = (2, 2, 2)
    population_list = [1, 1, 3]  # 1 Sr, 1 Ti, 3 O atoms
    species_list = [:Sr, :Ti, :O]
    charge_list = [+2, +4, -2]
    radii_list = [1.18, 0.42, 1.35]
    
    # Lattice parameters
    lattice_constant = 3.899  # Å
    L = lattice_constant * [1 0 0; 0 1 0; 0 0 1]
    lattice = Lattice(L, (true, true, true))
    
    # Ewald summation parameters
    depth = (4, 4, 4)
    alpha = 2 / lattice_constant
    
    # Build ion list and proximal pairs
    ion_list = build_ion_list(grid_size, species_list, charge_list, radii_list)
    proximal_pairs = build_proximal_pairs(ion_list, lattice, 0.75)
    
    # Solve the linear problem
    vector = build_vector(ion_list, lattice, interaction_energy, (alpha, depth, depth, depth))
    energy, solution_x, csp = build_linear_problem(grid_size, population_list, vector, proximal_pairs)
    
    selected_ions = []
    for index in CartesianIndices(ion_list)
        if solution_x[index] ≈ 1
            push!(selected_ions, ion_list[index])
        end
    end
    
    @test energy ≈ -6.061349350569214
    @test selected_ions == [Ion{3, Float64}(:Sr, 2, 1.18, [0.0, 0.5, 0.5]), Ion{3, Float64}(:Ti, 4, 0.42, [0.5, 0.0, 0.0]), Ion{3, Float64}(:O, -2, 1.35, [0.0, 0.0, 0.0]), Ion{3, Float64}(:O, -2, 1.35, [0.5, 0.5, 0.0]), Ion{3, Float64}(:O, -2, 1.35, [0.5, 0.0, 0.5])]

    # The quadratic problem formulation
    matrix = build_matrix(ion_list, lattice, interaction_energy, (alpha, depth, depth, depth))
    energy, solution_x, csp = build_quadratic_problem(grid_size, population_list, matrix)
    @test energy ≈ -6.061349350569214/2
end