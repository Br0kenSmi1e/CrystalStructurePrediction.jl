using CrystalStructurePrediction

"""
    setup_crystal_parameters()

Setup the parameters for SrTiO3 crystal structure prediction.
Returns grid_size, population_list, species_list, charge_list, radii_list, lattice, depth, and alpha.
"""
function setup_crystal_parameters()
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
    
    return grid_size, population_list, species_list, charge_list, radii_list, lattice, depth, alpha
end

"""
    run_crystal_structure_prediction(; use_quadratic_problem::Bool=false)

Run the crystal structure prediction for SrTiO3.

# Arguments
- `use_quadratic_problem::Bool`: Whether to use the quadratic problem. Only limited solvers support the quadratic problem, e.g. Gurobi (a commercial solver).
"""
function run_crystal_structure_prediction(; use_quadratic_problem::Bool=false)
    # Setup parameters
    grid_size, population_list, species_list, charge_list, radii_list, lattice, depth, alpha = setup_crystal_parameters()
    
    @info "Setting up crystal structure prediction for SrTiO3"
    @info "Grid size: $grid_size"
    @info "Population: $population_list $species_list"
    @info "Charges: $charge_list"
    @info "Ionic radii: $radii_list"
    
    # Build ion list and proximal pairs
    ion_list = build_ion_list(grid_size, species_list, charge_list, radii_list)
    @info "Created ion list with $(length(ion_list)) possible ion positions"
    
    proximal_pairs = build_proximal_pairs(ion_list, lattice, 0.75)
    @info "Identified $(length(proximal_pairs)) proximal pairs with cutoff 0.75"
    
    if use_quadratic_problem
        # Build interaction matrix
        @info "Building interaction energy matrix..."
        matrix = build_matrix(ion_list, lattice, interaction_energy, (alpha, depth, depth, depth))
    
        # Solve the quadratic problem
        @info "Solving quadratic optimization problem..."
        energy, solution_x, csp = build_quadratic_problem(grid_size, population_list, matrix, proximal_pairs)
    else
        # Build interaction vector
        @info "Building interaction energy vector..."
        vector = build_vector(ion_list, lattice, interaction_energy, (alpha, depth, depth, depth))
        
        # Solve the linear problem
        @info "Solving linear optimization problem..."
        energy, solution_x, csp = build_linear_problem(grid_size, population_list, vector, proximal_pairs)
    end
    
    # Display results
    @info "Optimization complete with energy: $energy"
    @info "Solution structure:"
    
    selected_ions = []
    for index in CartesianIndices(ion_list)
        if solution_x[index] ≈ 1
            push!(selected_ions, ion_list[index])
        end
    end
    
    for (i, ion) in enumerate(selected_ions)
        @info "Ion $i: $ion"
    end
    
    return energy, selected_ions, csp
end

# Run the prediction
energy, selected_ions, csp = run_crystal_structure_prediction()
