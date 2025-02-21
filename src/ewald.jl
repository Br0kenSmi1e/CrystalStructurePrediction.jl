# a great short note:
# https://www.cs.cornell.edu/courses/cs428/2006fa/Ewald%20Sum.pdf

# ewald sum = real space sum + reciprocal space sum + "self" energy

"""
    real_space_sum(depth, frac_pos_a, frac_pos_b, lattice, alpha)
Compute the real space part of Ewald summation.

Argument:
    lattice::AbstractMatrix: in the form [a, b, c], where a, b, c are column vectors.
"""
function real_space_sum(
        depth::Int,
        frac_pos_a::AbstractVector{T},
        frac_pos_b::AbstractVector{T},
        lattice::AbstractMatrix{T},
        alpha::T
        ) where T<:Real
    real_space_energy = zero(T)
    r_ab = lattice * (frac_pos_a - frac_pos_b)
    energy_in_the_cell = norm(r_ab) ≈ 0 ? 0 : erfc(sqrt(alpha) * norm(r_ab)) / (2 * norm(r_ab))
    real_space_energy += energy_in_the_cell
    for x in range(-depth, depth)
    for y in range(-depth, depth)
    for z in range(-depth, depth)
        if !(x==0 && y==0 && z==0)
            r = norm(r_ab + x * lattice[:, 1] + y * lattice[:, 2] + z * lattice[:, 3])
            real_space_energy += erfc(sqrt(alpha) * r) / (2 * r)
        end
    end
    end
    end
    return real_space_energy
end

function reciprocal_space_sum(
        depth::Int,
        frac_pos_a::AbstractVector{T},
        frac_pos_b::AbstractVector{T},
        lattice::AbstractMatrix{T},
        alpha::T
        ) where T
    
end