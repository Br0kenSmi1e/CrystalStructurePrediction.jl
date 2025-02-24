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