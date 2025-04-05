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

# (-6.061349350569213, Any[Ion{3, Float64}(:Sr, 2, 1.18, [0.5, 0.5, 0.0]), Ion{3, Float64}(:Ti, 4, 0.42, [0.0, 0.0, 0.5]), Ion{3, Float64}(:O, -2, 1.35, [0.0, 0.0, 0.0]), Ion{3, Float64}(:O, -2, 1.35, [0.5, 0.0, 0.5]), Ion{3, Float64}(:O, -2, 1.35, [0.0, 0.5, 0.5])], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0  …  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

# Visualize the crystal structure
using CairoMakie

function visualize_crystal_structure(selected_ions, lattice, shift)
    fig = Figure(; size = (300, 250))
    ax = Axis3(fig[1, 1], 
               aspect = :data,
               xlabel = "x", ylabel = "y", zlabel = "z",
               title = "SrTiO3 Crystal Structure")
       
    # Add unit cell edges
    cell_vertices = [
        [0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
        [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]
    ]
    
    # Convert to Cartesian coordinates
    cell_vertices_cart = [lattice.vectors * v for v in cell_vertices]
    
    # Define edges of the unit cell
    edges = [
        (1, 2), (2, 3), (3, 4), (4, 1),  # Bottom face
        (5, 6), (6, 7), (7, 8), (8, 5),  # Top face
        (1, 5), (2, 6), (3, 7), (4, 8)   # Connecting edges
    ]
    
    # Plot unit cell edges
    for (i, j) in edges
        lines!(ax, 
               [cell_vertices_cart[i][1], cell_vertices_cart[j][1]],
               [cell_vertices_cart[i][2], cell_vertices_cart[j][2]],
               [cell_vertices_cart[i][3], cell_vertices_cart[j][3]],
               color = :black, linewidth = 1)
    end
    
    # Plot the ions
    properties = Dict(:Sr => (color = :green, size = 30), :Ti => (color = :aqua, size = 25), :O => (color = :red, size = 10))
    
    # Plot each ion
    for ion in selected_ions
        # Plot the ion at its position and all periodic images within the unit cell
        for dx in -1:1, dy in -1:1, dz in -1:1
            # Add periodic image shift vector
            offset = [dx, dy, dz] .+ shift
            # Skip if the shifted position is outside the unit cell (0-1 range)
            shifted_pos = ion.frac_pos + offset
            if all(0 .<= shifted_pos .<= 1)
                # Convert to Cartesian coordinates
                shifted_cart_pos = lattice.vectors * shifted_pos
                scatter!(ax, [shifted_cart_pos[1]], [shifted_cart_pos[2]], [shifted_cart_pos[3]], 
                         color = properties[ion.species].color, 
                         markersize = properties[ion.species].size,
                         label = string(ion.species))
            end
        end
    end
 
    # Add legend with unique entries
    unique_species = unique([ion.species for ion in selected_ions])
    legend_elements = [MarkerElement(color = properties[sp].color, marker = :circle, markersize = properties[sp].size) for sp in unique_species]
    legend_labels = [string(sp) for sp in unique_species]
    
    Legend(fig[1, 2], legend_elements, legend_labels, "Species", patchsize = (30, 30))
    # Remove decorations and axis
    hidedecorations!(ax)
    hidespines!(ax)
    return fig
end

# Generate and save the visualization
lattice = setup_crystal_parameters()[6]
fig = visualize_crystal_structure(selected_ions, lattice, [0.5, 0.5, 0])

filename = joinpath(@__DIR__, "SrTiO3-structure.png")
save(filename, fig, dpi=20)
@info "Saved crystal structure visualization to: $filename"
