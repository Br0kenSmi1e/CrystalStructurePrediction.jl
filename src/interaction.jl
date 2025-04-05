# a great short note:
# https://www.cs.cornell.edu/courses/cs428/2006fa/Ewald%20Sum.pdf

# ewald sum = real space sum + reciprocal space sum + "self" energy
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

function real_space_potential(r::T, alpha::T) where T<:Real
    return r ≈ 0 ? - alpha / (π)^0.5 : erfc(alpha * r) / (2 * r)
end

"""
    real_space_sum(depth, frac_pos_a, frac_pos_b, lattice, alpha)
Compute the real space part of Ewald summation.

# Arguments
- `ion_a::Ion{D, T}`: the first ion.
- `ion_b::Ion{D, T}`: the second ion.
- `lattice::Lattice{D, T}`: the lattice.
- `alpha::T`: the Ewald parameter.
- `depth::NTuple{D, Int}`: summation depth on the dimensions.
"""
function real_space_sum(
        ion_a::Ion{D, T},
        ion_b::Ion{D, T},
        lattice::Lattice{D, T},
        alpha::T,
        depth::NTuple{D, Int}
        ) where {D, T}
    interaction = shift -> real_space_potential(norm(lattice.vectors * (ion_b.frac_pos + SVector(shift) - ion_a.frac_pos)), alpha)
    return charge(ion_a) * charge(ion_b) * 14.399645351950543 * periodic_sum(interaction, depth)
end

function reciprocal_space_potential(k::AbstractVector{T}, x::AbstractVector{T}, alpha::T) where T<:Real
    return norm(k) ≈ 0 ? 0 : (2π / dot(k, k)) * exp(-( dot(k, k) / (4*alpha^2) ) ) * cos(dot(k, x))
end

function reciprocal_space_sum(
        ion_a::Ion{D, T},
        ion_b::Ion{D, T},
        lattice::Lattice{D, T},
        alpha::T,
        depth::NTuple{D, Int}
        ) where {D, T}
    interaction = shift -> reciprocal_space_potential(2π * transpose(inv(lattice.vectors)) * SVector(shift), lattice.vectors * (ion_b.frac_pos - ion_a.frac_pos), alpha)
    return charge(ion_a) * charge(ion_b) * 14.399645351950543 * periodic_sum(interaction, depth) / abs(det(lattice.vectors))
end

# the minimum distance between two ions in a lattice
# TODO: this is the notoriously hard closest vector problem, try solve it with maybe integer programming?
function minimum_distance(pos_a::AbstractVector{T}, pos_b::AbstractVector{T}, lattice::Lattice{D, T}) where {D, T}
    return minimum(shift -> norm(lattice.vectors * (pos_b + SVector(shift) - pos_a)), Iterators.product(ntuple(x->-1:1, D)...))
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
    interaction = shift -> buckingham_potential(norm(lattice.vectors * (ion_b.frac_pos + SVector(shift) - ion_a.frac_pos)), buckingham_parameters(ion_a.type, ion_b.type)...)
    if ion_a ≈ ion_b
        return periodic_sum(interaction, depth)
    elseif minimum_distance(ion_a.frac_pos, ion_b.frac_pos, lattice) > threshold * (radii(ion_a) + radii(ion_b))
        return periodic_sum(interaction, depth)
    else
        return penalty
    end
end