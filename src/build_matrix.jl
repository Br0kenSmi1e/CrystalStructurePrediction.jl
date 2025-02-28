# function build_ion_sheet(grid_size::AbstractVector{Int},
#         species_list::AbstractVector{Symbol},
#         charge_list::AbstractVector{Int},
#         radii_list::AbstractVector{Float64}
#         )
#     ion_sheet = Array{Ion}(undef, length(species_list), prod(grid_size))
#     for (t, species) in enumerate(species_list)
#         for (p, pos) in enumerate(build_grid(grid_size))
#             ion_sheet[t, p] = Ion(species, charge_list[t], radii_list[t], pos ./ grid_size)
#         end
#     end
#     return ion_sheet
# end

function build_ion_list(grid_size::AbstractVector{Int},
        species_list::AbstractVector{Symbol},
        charge_list::AbstractVector{Int},
        radii_list::AbstractVector{Float64}
        )
    """
    return a list of ions in the format of:
    `[ion(t1, p1), ion(t1, p2), ..., ion(t1,pn), ion(t2, p1), ..., ion(tm, pn)]`
    """
    return [Ion(species_list[t], charge_list[t], radii_list[t], pos ./ grid_size) for t in range(1, length(species_list)) for pos in build_grid(grid_size)]
end

function interaction_energy(
        ion_a::Ion{T}, ion_b::Ion{T}, lattice::Lattice{T},
        alpha::T, real_depth::AbstractVector{Int}, reciprocal_depth::AbstractVector{Int}, buckingham_depth::AbstractVector{Int}
        ) where T<:Real
    return real_space_sum(ion_a, ion_b, lattice, alpha, real_depth) + reciprocal_space_sum(ion_a, ion_b, lattice, alpha, reciprocal_depth) + buckingham_sum(ion_a, ion_b, lattice, buckingham_depth)# + CrystalStructurePrediction.radii_penalty(ion_a, ion_b, lattice, 0.7)
end

function build_matrix(
        ion_list::AbstractVector{Ion{T}},
        lattice::Lattice{T},
        interaction_energy::FT,
        parameters::Tuple
        ) where {T<:Real, FT<:Function}
    matrix = zeros(T, length(ion_list), length(ion_list))
    for (index_a, ion_a) in enumerate(ion_list)
        for (index_b, ion_b) in enumerate(ion_list[index_a+1:end])
            matrix[index_a, index_a + index_b] = interaction_energy(ion_a, ion_b, lattice, parameters...)
        end
    end
    return (matrix + transpose(matrix)) / 2
end

function build_vector(
        ion_list::AbstractVector{Ion{T}},
        lattice::Lattice{T},
        interaction_energy::FT,
        parameters::Tuple
        ) where {T<:Real, FT<:Function}
    vector = zeros(T, length(ion_list) * (length(ion_list)-1) ÷ 2)
    for (index_a, ion_a) in enumerate(ion_list)
        for (index_b, ion_b) in enumerate(ion_list[index_a+1:end])
            vector[index_a + (index_a+index_b-1)*(index_a+index_b-2)÷2] = interaction_energy(ion_a, ion_b, lattice, parameters...)
        end
    end
    return vector
end

function build_proximal_pairs(ion_list::AbstractVector{Ion{T}}, lattice::Lattice{T}, c::Float64) where T<:Real
    isProximal = (i, j) -> CrystalStructurePrediction.minimum_distance(ion_list[i], ion_list[j], lattice) < c * (ion_list[i].radii + ion_list[j].radii)
    return [(i,j) for i in range(1, length(ion_list)) for j in range(i+1, length(ion_list)) if isProximal(i, j)]
end