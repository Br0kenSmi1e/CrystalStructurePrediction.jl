function build_ion_sheet(grid_size::AbstractVector{Int}, species_list::AbstractVector{Symbol}, charge_list::AbstractVector{Int})
    ion_sheet = Array{Ion}(undef, length(species_list), prod(grid_size))
    for (t, species) in enumerate(species_list)
        for (p, pos) in enumerate(build_grid(grid_size))
            ion_sheet[t, p] = Ion(species, charge_list[t], pos ./ grid_size)
        end
    end
    return ion_sheet
end

function interaction_energy(
        ion_a::Ion{T}, ion_b::Ion{T}, lattice::Lattice{T},
        alpha::T, real_depth::AbstractVector{Int}, reciprocal_depth::AbstractVector{Int}, buckingham_depth::AbstractVector{Int}
        ) where T<:Real
    return real_space_sum(ion_a, ion_b, lattice, alpha, real_depth) + reciprocal_space_sum(ion_a, ion_b, lattice, alpha, reciprocal_depth) + buckingham_sum(ion_a, ion_b, lattice, buckingham_depth)
end

function build_matrix(
        ion_sheet::AbstractMatrix{Ion},
        lattice::Lattice{T},
        interaction_energy::FT,
        parameters::Tuple
        ) where {T<:Real, FT<:Function}
    interaction_matrix = zeros(T, size(ion_sheet)..., size(ion_sheet)...)
    iter = CartesianIndices(ion_sheet)
    for (n, row_index) in enumerate(iter)
        for col_index in iter[n+1:end]
            interaction_matrix[row_index.I..., col_index.I...] = interaction_energy(ion_sheet[row_index], ion_sheet[col_index], lattice, parameters...)
        end
    end
    return interaction_matrix
end