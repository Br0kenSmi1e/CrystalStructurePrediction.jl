function build_ion_list(grid_size::AbstractVector{Int}, species_list::AbstractVector{String}, charge_list::AbstractVector{Int})
    return ion_list = [Ion(species, charge, pos ./ grid_size) for (species, charge) in zip(species_list, charge_list) for pos in build_grid(grid_size)]
end

function interaction_energy(
        ion_a::Ion{T}, ion_b::Ion{T}, lattice::Lattice{T},
        alpha::T, real_depth::AbstractVector{Int}, reciprocal_depth::AbstractVector{Int}, buckingham_depth::AbstractVector{Int}
        ) where T<:Real
    return real_space_sum(real_depth, ion_a, ion_b, lattice, alpha) + reciprocal_space_sum(reciprocal_depth, ion_a, ion_b, lattice, alpha) + buckingham_sum(buckingham_depth, ion_a, ion_b, lattice)
end

function build_matrix(
        ion_list::AbstractVector{Ion{T}},
        lattice::Lattice{T},
        alpha::T, real_depth::AbstractVector{Int}, reciprocal_depth::AbstractVector{Int}, buckingham_depth::AbstractVector{Int}
        ) where T<:Real
    interaction_matrix = Matrix{T}(undef, length(ion_list), length(ion_list))
    for (row, ion_a) in enumerate(ion_list)
        interaction_matrix[row, row] = zero(T)
        for (col, ion_b) in enumerate(ion_list[row+1:end])
            interaction_matrix[row, row + col] = interaction_energy(ion_a, ion_b, lattice, alpha, real_depth, reciprocal_depth, buckingham_depth)
            interaction_matrix[row + col, row] = interaction_matrix[row, row + col]
        end
    end
    return interaction_matrix
end