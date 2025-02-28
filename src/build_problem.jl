function build_problem(
        grid_size::AbstractVector{Int},
        population_list::AbstractVector{Int},
        interaction_matrix::AbstractMatrix{T}
        ) where T<:Real
    csp = Model(SCIP.Optimizer)
    num_grid_points = prod(grid_size)
    num_species = length(population_list)
    @variable(csp, x[1:num_species, 1:num_grid_points], Bin)
    @show interaction_matrix[48, 145]
    @objective(csp, Min, sum(interaction_matrix[i,j] * x[(i-1)÷num_grid_points + 1, (i-1)%num_grid_points + 1] * x[(j-1)÷num_grid_points + 1, (j-1)%num_grid_points + 1] for i in range(1, num_grid_points*num_species) for j in range(1, num_grid_points*num_species)))
    for i in range(1, num_species)
        @constraint(csp, sum(x[i, j] for j in range(1, num_grid_points)) == population_list[i])
    end
    for j in range(1, num_grid_points)
        @constraint(csp, sum(x[i, j] for i in range(1, num_species)) <= 1)
    end
    optimize!(csp)
    assert_is_solved_and_feasible(csp)
    return value.(x)
end

function build_linear_problem(
        grid_size::AbstractVector{Int},
        population_list::AbstractVector{Int},
        interaction_vector::AbstractVector{T},
        proximal_pairs::AbstractVector{Tuple{Int, Int}}
        ) where T<:Real
    csp = Model(Gurobi.Optimizer)
    num_grid_points = prod(grid_size)
    num_species = length(population_list)
    @variable(csp, 0 <= x[1:num_species*num_grid_points] <= 1, Int)
    @variable(csp, 0 <= s[1:num_species*num_grid_points*(num_species*num_grid_points-1)÷2] <= 1, Int)
    # @constraint(csp, sum(x[t, 1] for t in range(1, num_species)) == 1)
    for t in range(1, num_species)
        @constraint(csp, sum(x[num_grid_points*(t-1)+p] for p in range(1, num_grid_points)) == population_list[t])
    end
    for p in range(1, num_grid_points)
        @constraint(csp, sum(x[num_grid_points*(t-1)+p] for t in range(1, num_species)) <= 1)
    end
    for (i, j) in proximal_pairs
        @constraint(csp, x[i] + x[j] <= 1)
    end
    for i in range(1, num_species*num_grid_points-1)
    for j in range(i+1, num_species*num_grid_points)
        @constraint(csp, s[i + (j-1)*(j-2)÷2] <= x[i])
        @constraint(csp, s[i + (j-1)*(j-2)÷2] <= x[j])
        @constraint(csp, s[i + (j-1)*(j-2)÷2] >= x[i] + x[j] - 1)
    end
    end 
    @objective(csp, Min, dot(interaction_vector, s))
    optimize!(csp)
    assert_is_solved_and_feasible(csp)
    return objective_value(csp), value.(x), value.(s)
end

function build_quadratic_problem(
        grid_size::AbstractVector{Int},
        population_list::AbstractVector{Int},
        interaction_matrix::AbstractMatrix{T},
        proximal_pairs::AbstractVector{Tuple{Int, Int}}
        ) where T<:Real

    csp = Model(Gurobi.Optimizer)
    num_grid_points = prod(grid_size)
    num_species = length(population_list)
    @variable(csp, 0 <= x[1:num_species*num_grid_points] <= 1, Int)
    # @constraint(csp, sum(x[t, 1] for t in range(1, num_species)) == 1)
    for t in range(1, num_species)
        @constraint(csp, sum(x[num_grid_points*(t-1)+p] for p in range(1, num_grid_points)) == population_list[t])
    end
    for p in range(1, num_grid_points)
        @constraint(csp, sum(x[num_grid_points*(t-1)+p] for t in range(1, num_species)) <= 1)
    end
    for (i, j) in proximal_pairs
        @constraint(csp, x[i] + x[j] <= 1)
    end
    # for index_a in CartesianIndices(ion_sheet)
    #     for index_b in CartesianIndices(ion_sheet)
    #         if norm(lattice.vectors * (ion_sheet[index_a].frac_pos - ion_sheet[index_b].frac_pos)) < 0.75 * (ion_sheet[index_a].radii + ion_sheet[index_b].radii)
    #             @constraint(csp, x[index_a] + x[index_b] <= 1)
    #         end
    #     end
    # end
    @objective(csp, Min, x' * interaction_matrix * x)
    optimize!(csp)
    assert_is_solved_and_feasible(csp)
    return objective_value(csp), value.(x)
end