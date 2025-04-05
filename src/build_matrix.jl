"""
return a list of ions in the format of:
`[ion(t1, p1), ion(t1, p2), ..., ion(t1,pn), ion(t2, p1), ..., ion(tm, pn)]`
"""
function build_ion_list(grid_size::NTuple{N, Int},
        type_list::AbstractVector{IonType{T}},
        ) where {N, T}
    return [Ion(type_list[t], (ci.I .- 1) ./ grid_size) for t in range(1, length(type_list)) for ci in vec(CartesianIndices(grid_size))]
end

function build_matrix(
        ion_list::AbstractVector{Ion{D, T}},
        lattice::Lattice{D, T},
        interaction_energy::FT,
        parameters::Tuple
        ) where {D, T<:Real, FT<:Function}
    matrix = zeros(T, length(ion_list), length(ion_list))
    for (index_a, ion_a) in enumerate(ion_list)
        for (index_b, ion_b) in enumerate(ion_list[index_a:end])
            matrix[index_a, index_a + index_b - 1] = interaction_energy(ion_a, ion_b, lattice, parameters...)
        end
    end
    return (matrix + transpose(matrix)) / 2
end

function build_vector(
        ion_list::AbstractVector{Ion{D, T}},
        lattice::Lattice{D, T},
        interaction_energy::FT,
        parameters::Tuple
        ) where {D, T<:Real, FT<:Function}
    vector = zeros(T, length(ion_list) * (length(ion_list)-1) ÷ 2)
    for (index_a, ion_a) in enumerate(ion_list)
        for (index_b, ion_b) in enumerate(ion_list[index_a+1:end])
            vector[index_a + (index_a+index_b-1)*(index_a+index_b-2)÷2] = interaction_energy(ion_a, ion_b, lattice, parameters...)
        end
    end
    return vector
end

function build_proximal_pairs(ion_list::AbstractVector{Ion{D, T}}, lattice::Lattice{D, T}, c::Float64) where {D, T}
    isProximal = (i, j) -> CrystalStructurePrediction.minimum_distance(ion_list[i].frac_pos, ion_list[j].frac_pos, lattice) < c * (radii(ion_list[i]) + radii(ion_list[j]))
    return [(i,j) for i in range(1, length(ion_list)) for j in range(i+1, length(ion_list)) if isProximal(i, j)]
end