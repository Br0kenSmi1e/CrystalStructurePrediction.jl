struct IonOptimizationResult{D, T}
    energy::T
    selected_ions::AbstractVector{Ion{D, T}}
end

"""
    optimize_linear(interaction, ions, populations, lattice; optimizer = SCIP.Optimizer, optimizer_options = Dict(), proximal_threshold = 0.75)

Build a linear problem for crystal structure prediction. Suited for solvers not supporting quadratic constraints.

# Arguments
- `interaction::Function`: The interaction function.
- `ions::AbstractVector{Ion{D, T}}`: The ions.
- `populations::Dict{IonType{T}, Int}`: The populations of the ions.
- `lattice::Lattice`: The lattice.

# Keyword arguments
- `optimizer`: The optimizer.
- `optimizer_options`: The options for the optimizer, e.g. `optimizer_options = Dict("NodefileSave" => 1)` for Gurobi.
- `proximal_threshold`: The threshold for the proximal pairs.
"""
function optimize_linear(interaction,
            ions::AbstractVector{Ion{D, T}},
            populations::Dict{IonType{T}, Int},  # a dictionary of ion types and their populations
            lattice::Lattice;
            optimizer = SCIP.Optimizer,
            optimizer_options = Dict(),
            proximal_threshold = 0.75
        ) where {D, T<:Real}
    # build the model
    csp = Model(optimizer); set_silent(csp)
    num_ions = length(ions)
    unique_coordinates = unique!([ion.frac_pos for ion in ions])
    @variable(csp, 0 <= x[1:num_ions] <= 1, Int)
    @variable(csp, 0 <= s[1:num_ions*(num_ions-1)÷2] <= 1, Int)

    # each ion type has a constraint on the number of ions of that type
    for (type, population) in populations
        @constraint(csp, sum(x[i] for i in range(1, num_ions) if ions[i].type == type) == population)
    end
    # each grid point has at most one ion
    for coord in unique_coordinates
        @constraint(csp, sum(x[i] for i in range(1, num_ions) if ions[i].frac_pos == coord) <= 1)
    end
    # ions are not allowed to be too close to each other
    for i in 1:num_ions-1, j in i+1:num_ions
        if too_close(ions[i], ions[j], lattice, proximal_threshold)
            @constraint(csp, x[i] + x[j] <= 1)
        end
    end
    # s is a "matrix" representing the co-existence of ions
    energy = zero(eltype(s))
    for i in range(1, num_ions-1)
        for j in range(i+1, num_ions)
            @constraint(csp, s[i + (j-1)*(j-2)÷2] <= x[i])
            @constraint(csp, s[i + (j-1)*(j-2)÷2] <= x[j])
            @constraint(csp, s[i + (j-1)*(j-2)÷2] >= x[i] + x[j] - 1)
            energy += interaction(ions[i], ions[j], lattice) * s[i + (j-1)*(j-2)÷2]
        end
    end 
    # minimize the interaction energy
    @objective(csp, Min, energy)
    # set the optimizer options
    for (key, value) in optimizer_options
        set_optimizer_attribute(csp, key, value)
    end
    # solve the problem
    optimize!(csp)
    assert_is_solved_and_feasible(csp)
    return IonOptimizationResult(objective_value(csp), [ions[i] for i in 1:num_ions if value.(x)[i] ≈ 1])
end

function too_close(ion_a::Ion, ion_b::Ion, lattice::Lattice, c::Float64)
    return minimum_distance(ion_a.frac_pos, ion_b.frac_pos, lattice) < c * (radii(ion_a) + radii(ion_b))
end

"""
    optimize_quadratic(interaction, ions, populations, lattice; optimizer = SCIP.Optimizer, optimizer_options = Dict(), proximal_threshold = 0.75)

Build a quadratic problem for crystal structure prediction.

# Arguments
- `interaction::Function`: The interaction function.
- `ions::AbstractVector{Ion{D, T}}`: The ions.
- `populations::Dict{IonType{T}, Int}`: The populations of the ions.
- `lattice::Lattice`: The lattice.

# Keyword arguments
- `optimizer`: The optimizer.
- `optimizer_options`: The options for the optimizer, e.g. `optimizer_options = Dict("NodefileSave" => 1)` for Gurobi.
- `proximal_threshold`: The threshold for the proximal pairs.
"""
function optimize_quadratic(interaction,
            ion_list::AbstractVector{Ion{D, T}},
            populations::Dict{IonType{T}, Int},
            lattice::Lattice;
            optimizer = SCIP.Optimizer,
            optimizer_options = Dict(),
            proximal_threshold = 0.75
        ) where {D, T<:Real}
    # build the model
    csp = Model(optimizer); set_silent(csp)
    num_ions = length(ion_list)
    unique_coordinates = unique!([ion.frac_pos for ion in ion_list])
    @variable(csp, x[1:num_ions], Bin)
    # each ion type has a constraint on the number of ions of that type
    for (type, population) in populations
        @constraint(csp, sum(x[i] for i in range(1, num_ions) if ion_list[i].type == type) == population)
    end
    # each grid point has at most one ion
    for coord in unique_coordinates
        @constraint(csp, sum(x[i] for i in range(1, num_ions) if ion_list[i].frac_pos == coord) <= 1)
    end
    # ions are not allowed to be too close to each other
    for i in 1:num_ions-1, j in i+1:num_ions
        if too_close(ion_list[i], ion_list[j], lattice, proximal_threshold)
            @constraint(csp, x[i] + x[j] <= 1)
        end
    end
    # minimize the interaction energy
    @objective(csp, Min, sum(interaction(ion_list[i], ion_list[j], lattice)*x[i]*x[j] for i in range(1, num_ions) for j in range(i+1, num_ions)))
    # set the optimizer options
    for (key, value) in optimizer_options
        set_optimizer_attribute(csp, key, value)
    end
    # solve the problem
    optimize!(csp)
    assert_is_solved_and_feasible(csp)
    return IonOptimizationResult(objective_value(csp), [ion_list[i] for i in 1:num_ions if value.(x)[i] ≈ 1])
end

