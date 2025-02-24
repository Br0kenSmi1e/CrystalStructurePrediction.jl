# a great short note:
# https://www.cs.cornell.edu/courses/cs428/2006fa/Ewald%20Sum.pdf

# ewald sum = real space sum + reciprocal space sum + "self" energy


function build_shifts(depth::AbstractVector{Int})
    return [shift - depth for shift in build_grid(2 * depth .+ 1)]
end

"""
    summation(depth, interaction)
Compute periodic summation of interaction to given depth.

Argument:
    depth::AbstractVector{Int}: summation depth on the dimensions.
    interaction<:Function: the interaction energy, given the shift as a list of integer.
"""
function summation(depth::AbstractVector{Int}, interaction::FT) where FT<:Function
    energy = 0
    for shift in build_shifts(depth)
        energy += interaction(shift)
    end
    return energy
end

function real_space_potential(r::T, alpha::T) where T<:Real
    return r ≈ 0 ? 0 : erfc(alpha * r) / (2 * r)
end

"""
    real_space_sum(depth, frac_pos_a, frac_pos_b, lattice, alpha)
Compute the real space part of Ewald summation.

Argument:
    lattice::AbstractMatrix: in the form [a, b, c], where a, b, c are column vectors.
"""
function real_space_sum(
        depth::AbstractVector{Int},
        ion_a::Ion{T},
        ion_b::Ion{T},
        lattice::Lattice{T},
        alpha::T
        ) where T<:Real
    interaction = shift -> real_space_potential(norm(lattice.vectors * (ion_b.frac_pos + shift - ion_a.frac_pos)), alpha)
    return ion_a.charge * ion_b.charge * summation(depth, interaction)
end

# function reciprocal_space_potential(
#         shift::AbstractVector{Int},
#         frac_field::AbstractVector{T},
#         frac_source::AbstractVector{T},
#         lattice::Lattice{T},
#         alpha::T
#         ) where T<:Real
#     k = 2π * transpose(inv(lattice.vectors)) * shift
#     return norm(k) ≈ 0 ? 0 : (2π / norm(k)^2) * exp(-( norm(k)^2 / (4*alpha^2) ) + 2π*im* dot(shift, frac_source - frac_field))
# end

function reciprocal_space_potential(k::AbstractVector{T}, x::AbstractVector{T}, alpha::T) where T<:Real
    return norm(k) ≈ 0 ? 0 : (2π / norm(k)^2) * exp(-( norm(k)^2 / (4*alpha^2) ) + 2π*im* dot(k, x))
end

function reciprocal_space_sum(
        depth::AbstractVector{Int},
        ion_a::Ion{T},
        ion_b::Ion{T},
        lattice::Lattice{T},
        alpha::T
        ) where T<:Real
    interaction = shift -> reciprocal_space_potential(2π * transpose(inv(lattice.vectors)) * shift, lattice.vectors * (ion_b.frac_pos - ion_a.frac_pos), alpha)
    return ion_a.charge * ion_b.charge * summation(depth, interaction) / abs(det(lattice.vectors))
end

function buckingham_potential(r::T, A::T, ρ::T, C::T) where T<:Real
    return A * exp(-r/ρ) - C / r^6
end

function buckingham_parameters(ion_a::Ion{T}, ion_b::Ion{T}) where T
    if ion_a.species=="O" && ion_b.species=="O"
        return 1388.7, 0.36262, 175.0
    elseif ion_a.species=="Sr" && ion_b.species=="O"
        return 1952.39, 0.33685, 19.22
    elseif ion_a.species=="Ti" && ion_b.species=="O"
        return 4590.7279, 0.261, 0.0
    elseif ion_a.species=="O" && ion_b.species=="Sr"
        return 1952.39, 0.33685, 19.22
    elseif ion_a.species=="O" && ion_b.species=="Ti"
        return 4590.7279, 0.261, 0.0
    else
        return 0.0, 1.0, 0.0
    end
end

function buckingham_sum(
        depth::AbstractVector{Int},
        ion_a::Ion{T},
        ion_b::Ion{T},
        lattice::Lattice{T}
        ) where T<:Real
    interaction = shift -> buckingham_potential(norm(lattice.vectors * (ion_b.frac_pos + shift - ion_a.frac_pos)), buckingham_parameters(ion_a, ion_b)...)
    return summation(depth, interaction)
end