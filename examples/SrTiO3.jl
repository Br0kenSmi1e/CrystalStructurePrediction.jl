using CrystalStructurePrediction
using CairoMakie, SCIP

# Run the crystal structure prediction, alpha is the Ewald parameter
function run_crystal_structure_prediction(grid_size, populations, lattice; use_quadratic_problem::Bool=false)
    alpha = 2 / maximum(lattice.vectors)
    @info "Setting up crystal structure prediction with"
    @info "Grid size: $grid_size"
    @info "Populations: $populations"
    @info "Ewald parameter: $alpha"
    
    # Build ion list and proximal pairs
    ion_list = ions_on_grid(grid_size, collect(keys(populations)))
    @info "Created ion list with $(length(ion_list)) possible ion positions"
    
    # Ewald summation parameters
    depth = (4, 4, 4)
    if use_quadratic_problem
        # Solve with the quadratic formulation
        @info "Solving quadratic optimization problem..."
        res = optimize_quadratic(ion_list, populations, lattice; optimizer=SCIP.Optimizer) do ion_a, ion_b, lattice
            interaction_energy(ion_a, ion_b, lattice, alpha, depth, depth, depth)
        end
    else
        # Solve with the linear formulation
        @info "Solving linear optimization problem..."
        res = optimize_linear(ion_list, populations, lattice; optimizer=SCIP.Optimizer) do ion_a, ion_b, lattice
            interaction_energy(ion_a, ion_b, lattice, alpha, depth, depth, depth)
        end
    end
    # Display results
    @info "Optimization complete with energy: $(res.energy)"
    for ion in res.selected_ions
        @info "Ion: $ion"
    end
    return res
end

# Visualize the crystal structure
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
    properties = Dict(:Sr => (color = :green,), :Ti => (color = :aqua,), :O => (color = :red,))
    
    # Plot each ion
    for ion in selected_ions
        # Plot the ion at its position and all periodic images within the unit cell
        coordinates = []
        for dx in -1:1, dy in -1:1, dz in -1:1
            # Add periodic image shift vector
            offset = [dx, dy, dz] .+ shift
            # Skip if the shifted position is outside the unit cell (0-1 range)
            shifted_pos = ion.frac_pos + offset
            if all(0 .<= shifted_pos .<= 1)
                # Convert to Cartesian coordinates
                shifted_cart_pos = lattice.vectors * shifted_pos
                push!(coordinates, shifted_cart_pos)
            end
        end
        scatter!(ax, coordinates, 
                    color = properties[ion.type.species].color, 
                    markersize = ion.type.radii * 20,
                    label = string(ion.type.species))
    end
 
    # Add legend with unique entries
    unique_species = unique([ion.type for ion in selected_ions])
    legend_elements = [MarkerElement(color = properties[sp.species].color, marker = :circle, markersize = sp.radii * 20) for sp in unique_species]
    legend_labels = [string(sp.species) for sp in unique_species]
    
    Legend(fig[1, 2], legend_elements, legend_labels, "Species", patchsize = (30, 30))
    # Remove decorations and axis
    hidedecorations!(ax)
    hidespines!(ax)
    return fig
end

####### Run the prediction #######
function run_SrTiO3_prediction()
    # Crystal structure parameters
    lattice_constant = 3.899  # Å
    lattice = Lattice(lattice_constant .* [1 0 0; 0 1 0; 0 0 1])
    grid_size = (2, 2, 2)
    populations = Dict(
        IonType(:Sr, +2, 1.18) => 1,  # 1 Sr atom
        IonType(:Ti, +4, 0.42) => 1,  # 1 Ti atom
        IonType(:O, -2, 1.35) => 3    # 3 O atoms
    )

    res = run_crystal_structure_prediction(grid_size, populations, lattice; use_quadratic_problem=false)

    # Generate and save the visualization
    origin = res.selected_ions[findfirst(x -> x.type.species == :Sr, res.selected_ions)].frac_pos
    fig = visualize_crystal_structure(res.selected_ions, lattice, origin)

    filename = joinpath(@__DIR__, "SrTiO3-structure.png")
    save(filename, fig, dpi=20)
    @info "Saved crystal structure visualization to: $filename"
end

run_SrTiO3_prediction()
