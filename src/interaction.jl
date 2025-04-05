# a great short note:
# https://www.cs.cornell.edu/courses/cs428/2006fa/Ewald%20Sum.pdf

# ewald sum = real space sum + reciprocal space sum + "self" energy

"""
    real_space_Ewald_sum(ion_a, ion_b, lattice, alpha, depth)

Compute the real space part of Ewald summation.

# Arguments
- `ion_a::Ion{D, T}`: the first ion.
- `ion_b::Ion{D, T}`: the second ion.
- `lattice::Lattice{D, T}`: the lattice.
- `alpha::T`: the Ewald parameter.
- `depth::NTuple{D, Int}`: summation depth on the dimensions.
"""
function interaction_energy(
        ion_a::Ion{D, T}, ion_b::Ion{D, T}, lattice::Lattice{D, T},
        alpha::T, real_depth::NTuple{D, Int}, reciprocal_depth::NTuple{D, Int}, buckingham_depth::NTuple{D, Int}
        ) where {D, T}
    return real_space_Ewald_sum(ion_a, ion_b, lattice, alpha, real_depth) + reciprocal_space_Ewald_sum(ion_a, ion_b, lattice, alpha, reciprocal_depth) + buckingham_sum(ion_a, ion_b, lattice, buckingham_depth)# + radii_penalty(ion_a, ion_b, lattice, 0.75)

    # Q: What is the constant 14.399645351950543?
    real_space_energy = charge(ion_a) * charge(ion_b) * 14.399645351950543 * periodic_sum(real_depth) do shift
        real_space_potential(distance(lattice, ion_b.frac_pos + SVector(shift), ion_a.frac_pos), alpha)
    end
    reciprocal_space_energy = charge(ion_a) * charge(ion_b) * 14.399645351950543 * periodic_sum(reciprocal_depth) do shift
        reciprocal_space_potential(reciprocal_vectors(lattice) * SVector(shift), cartesian(lattice, ion_b.frac_pos - ion_a.frac_pos), alpha)
    end
    buckingham_energy = charge(ion_a) * charge(ion_b) * 14.399645351950543 * periodic_sum(buckingham_depth) do shift
        buckingham_potential(norm(cartesian(lattice, ion_b.frac_pos + SVector(shift) - ion_a.frac_pos)), buckingham_parameters(ion_a.type, ion_b.type)...)
    end
    return real_space_energy + reciprocal_space_energy + buckingham_energy
end

function real_space_Ewald_sum(ion_a::Ion{D, T}, ion_b::Ion{D, T}, lattice::Lattice{D, T}, alpha::T, depth::NTuple{D, Int}) where {D, T}
    interaction = shift -> real_space_potential(distance(lattice, ion_b.frac_pos + SVector(shift), ion_a.frac_pos), alpha)
    return charge(ion_a) * charge(ion_b) * 14.399645351950543 * periodic_sum(interaction, depth)
end

# What is this potential function? Does it have a name?
function real_space_potential(r::T, alpha::T) where T<:Real
    return T(r ≈ 0 ? - alpha / sqrt(π) : erfc(alpha * r) / (2 * r))
end

function reciprocal_space_potential(k::AbstractVector{T}, x::AbstractVector{T}, alpha::T) where T<:Real
    return norm(k) ≈ 0 ? 0 : (2π / dot(k, k)) * exp(-( dot(k, k) / (4*alpha^2) ) ) * cos(dot(k, x))
end

function reciprocal_space_Ewald_sum(
        ion_a::Ion{D, T},
        ion_b::Ion{D, T},
        lattice::Lattice{D, T},
        alpha::T,
        depth::NTuple{D, Int}
        ) where {D, T}
    interaction = shift -> reciprocal_space_potential(reciprocal_vectors(lattice) * SVector(shift), cartesian(lattice, ion_b.frac_pos - ion_a.frac_pos), alpha)
    return charge(ion_a) * charge(ion_b) * 14.399645351950543 * periodic_sum(interaction, depth) / abs(det(lattice.vectors))
end

"""
    periodic_sum(interaction, depth)

Compute periodic summation of interaction to given depth.

# Arguments
- `interaction<:Function`: the interaction energy, given the shift as a list of integers.
- `depth::NTuple{N, Int}`: summation depth on the dimensions.
"""
function periodic_sum(interaction::FT, depth::NTuple{N, Int}) where {N, FT<:Function}
    return sum(interaction, Iterators.product(ntuple(i->-depth[i]:depth[i], N)...))
end

function radii_penalty(ion_a::Ion{D, T}, ion_b::Ion{D, T}, lattice::Lattice{D, T}, c::T, penalty::T=3e2) where {D, T}
    if ion_a ≈ ion_b
        return zero(T)
    else
        return minimum_distance(ion_a.frac_pos, ion_b.frac_pos, lattice) / (radii(ion_a) + radii(ion_b)) > c ? zero(T) : penalty
    end
end

function buckingham_potential(r::T, A::T, ρ::T, C::T) where T<:Real
    return r ≈ 0 ? 0 : A * exp(-r/ρ) - C / r^6
end

function buckingham_parameters(ion_a::IonType{T}, ion_b::IonType{T}) where T<:Real
    if Set([ion_a.species, ion_b.species]) == Set([:O])
        return 1388.7, 0.36262, 175.0
    elseif Set([ion_a.species, ion_b.species]) == Set([:Sr, :O])
        return 1952.39, 0.33685, 19.22
    elseif Set([ion_a.species, ion_b.species]) == Set([:Ti, :O])
        return 4590.7279, 0.261, 0.0
    else
        # TODO: maybe throw an error here?
        return 0.0, 1.0, 0.0
    end
end

function buckingham_sum(
        ion_a::Ion{D, T},
        ion_b::Ion{D, T},
        lattice::Lattice{D, T},
        depth::NTuple{D, Int},
        threshold::Float64=0.75,
        penalty::Float64=3e2
        ) where {D, T}
    interaction = shift -> buckingham_potential(norm(cartesian(lattice, ion_b.frac_pos + SVector(shift) - ion_a.frac_pos)), buckingham_parameters(ion_a.type, ion_b.type)...)
    if ion_a ≈ ion_b
        return periodic_sum(interaction, depth)
    elseif minimum_distance(ion_a.frac_pos, ion_b.frac_pos, lattice) > threshold * (radii(ion_a) + radii(ion_b))
        return periodic_sum(interaction, depth)
    else
        return penalty
    end
end