"""
    BUCKINGHAM_PARAMETERS

A map from the species of two ions to the Buckingham parameters: `(A, ρ, C)`.
- `A`: the well depth.
- `ρ`: the range parameter.
- `C`: the Buckingham constant.
"""
const BUCKINGHAM_PARAMETERS = Dict(
    minmax(:O, :O) => (1388.7, 0.36262, 175.0),
    minmax(:Sr, :O) => (1952.39, 0.33685, 19.22),
    minmax(:Ti, :O) => (4590.7279, 0.261, 0.0)
)

"""
    interaction_energy(ion_a, ion_b, lattice, alpha, real_depth, reciprocal_depth, buckingham_depth, buckingham_threshold, buckingham_penalty, radii_threshold, radii_penalty)

Compute the interaction energy between two ions, which is a sum of the real space Ewald sum, the reciprocal space Ewald sum, the Buckingham potential, and the radii penalty.

# Arguments
- `ion_a::Ion{D, T}`: the first ion.
- `ion_b::Ion{D, T}`: the second ion.
- `lattice::Lattice{D, T}`: the lattice.
- `alpha::T`: the Ewald parameter.
- `real_depth::NTuple{D, Int}`: summation depth on the dimensions for the real space Ewald sum.
- `reciprocal_depth::NTuple{D, Int}`: summation depth on the dimensions for the reciprocal space Ewald sum.
- `buckingham_depth::NTuple{D, Int}`: summation depth on the dimensions for the Buckingham potential.

# Keyword arguments
- `buckingham_threshold::T=0.75`: the threshold for the Buckingham potential.
- `buckingham_penalty::T=3e2`: the penalty for the Buckingham potential.
- `radii_threshold::T=0.0`: the threshold for the radii penalty, defined as the ratio of the distance between the two ions to the sum of their radii.
- `radii_penalty::T=0.0`: the penalty for the radii penalty.

# Notes
- The Buckingham potential is only summed over the species pairs that are present in the `CrystalStructurePrediction.BUCKINGHAM_PARAMETERS` dictionary:
  $BUCKINGHAM_PARAMETERS
  If you want to use a different potential, you can do so by setting the `BUCKINGHAM_PARAMETERS` dictionary to your desired potential.

# References
- a great short note: https://www.cs.cornell.edu/courses/cs428/2006fa/Ewald%20Sum.pdf
"""
function interaction_energy(
        ion_a::Ion{D, T}, ion_b::Ion{D, T}, lattice::Lattice{D, T},
        alpha::T, real_depth::NTuple{D, Int}, reciprocal_depth::NTuple{D, Int}, buckingham_depth::NTuple{D, Int},
        buckingham_threshold::T=0.75, buckingham_penalty::T=3e2, radii_threshold::T=0.0, radii_penalty::T=0.0
        ) where {D, T}
    # k * e / Ang^2 = 14.399645351950543, where k is the Coulomb constant, e is the elementary charge, and Ang is the Angstrom unit 10^-10 m.
    # real space Ewald sum
    energy = zero(T)
    energy += charge(ion_a) * charge(ion_b) * 14.399645351950543 * periodic_sum(real_depth) do shift
        real_space_potential(distance(lattice, ion_b.frac_pos + SVector(shift), ion_a.frac_pos), alpha)
    end
    # reciprocal space Ewald sum
    energy += charge(ion_a) * charge(ion_b) * 14.399645351950543 / volume(lattice) * periodic_sum(reciprocal_depth) do shift
        k = reciprocal_vectors(lattice) * SVector(shift)
        norm(k) ≈ 0 && return zero(T)  # only sum over non-zero k
        reciprocal_space_potential(k, cartesian(lattice, ion_b.frac_pos - ion_a.frac_pos), alpha)
    end
    # buckingham energy
    buckingham_key = minmax(species(ion_a), species(ion_b))
    if haskey(BUCKINGHAM_PARAMETERS, buckingham_key)
        energy += if (ion_a == ion_b || minimum_distance(ion_a.frac_pos, ion_b.frac_pos, lattice) > buckingham_threshold * (radii(ion_a) + radii(ion_b)))
            params = BUCKINGHAM_PARAMETERS[buckingham_key]
            periodic_sum(buckingham_depth) do shift
                r = distance(lattice, ion_b.frac_pos + SVector(shift), ion_a.frac_pos)
                r ≈ 0 && return zero(T)  # only sum over non-zero r
                buckingham_potential(r, params...)
            end
        else
            buckingham_penalty
        end
    end
    # radii penalty
    if !(ion_a == ion_b) && minimum_distance(ion_a.frac_pos, ion_b.frac_pos, lattice) / (radii(ion_a) + radii(ion_b)) > radii_threshold
        energy += radii_penalty
    end
    return energy
end

# What is this potential function? Does it have a name?
function real_space_potential(r::T, alpha::T) where T<:Real
    return T(r ≈ 0 ? -alpha / sqrt(π) : erfc(alpha * r) / (2 * r))
end

function reciprocal_space_potential(k::AbstractVector{T}, x::AbstractVector{T}, alpha::T) where T<:Real
    return 2π / dot(k, k) * exp(-( dot(k, k) / (4*alpha^2) ) ) * cos(dot(k, x))
end

function buckingham_potential(r::T, A::T, ρ::T, C::T) where T<:Real
    return A * exp(-r/ρ) - C / r^6
end

"""
    ions_on_grid(grid_size::NTuple{N, Int}, type_list::AbstractVector{IonType{T}}) where {N, T}

Create a list of ions on a grid.

# Arguments
- `grid_size::NTuple{N, Int}`: the size of the grid.
- `type_list::AbstractVector{IonType{T}}`: the list of ion types.

# Returns
- A vector of ions: `[ion(t1, p1), ion(t1, p2), ..., ion(t1,pn), ion(t2, p1), ..., ion(tm, pn)]`.
where `t1, ..., tm` are the types of the ions and `p1, ..., pn` are the fractional positions of the ions on the grid.
"""
function ions_on_grid(grid_size::NTuple{N, Int}, type_list::AbstractVector{IonType{T}}) where {N, T}
    return [Ion(type_list[t], SVector((ci.I .- 1) .// grid_size)) for t in range(1, length(type_list)) for ci in CartesianIndices(grid_size)]
end