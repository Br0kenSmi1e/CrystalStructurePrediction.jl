# a great short note:
# https://www.cs.cornell.edu/courses/cs428/2006fa/Ewald%20Sum.pdf

# ewald sum = real space sum + reciprocal space sum + "self" energy

struct Lattice{T}
    vectors::AbstractMatrix{T}
    pbc::AbstractVector{Bool}
end

import Base.size
size(lattice::Lattice) = size(lattice.vectors)

function periodic_vectors(lattice::Lattice{T}) where T
    periodic_vectors = Matrix{T}(undef, size(lattice)[1], sum(lattice.pbc))
    for column in range(1, size(lattice)[2])
        if lattice.pbc[column]
            periodic_vectors[:, column] = lattice.vectors[:, column]
        end
    end
    return periodic_vectors
end

function build_shifts(depth::AbstractVector{Int})
    # shifts = Array{NTuple{length(depth), Int}}(undef, Tuple(2 * depth .+ 1))
    function decompose(n::Int)
        remain = n
        result = Vector{Int}(undef, length(depth))
        for i in range(1, length(depth))
            result[i] = remain ÷ prod(2 * depth[i+1:end] .+ 1) - depth[i]
            remain = remain % prod(2 * depth[i+1:end] .+ 1)
        end
        return result
    end
    return [decompose(n) for n in range(0, prod(2 * depth .+ 1) - 1) if n != (prod(2 * depth .+ 1) - 1) ÷ 2]
end

"""
    real_space_sum(depth, frac_pos_a, frac_pos_b, lattice, alpha)
Compute the real space part of Ewald summation.

Argument:
    lattice::AbstractMatrix: in the form [a, b, c], where a, b, c are column vectors.
"""
function real_space_sum(
        depth::AbstractVector{Int},
        frac_pos_a::AbstractVector{T},
        frac_pos_b::AbstractVector{T},
        lattice::Lattice{T},
        alpha::T
        ) where T<:Real
    r_ab = lattice.vectors * (frac_pos_a - frac_pos_b)
    real_space_energy = norm(r_ab) ≈ 0 ? 0 : erfc(alpha * norm(r_ab)) / (2 * norm(r_ab))
    for shift in build_shifts(depth)
        r = norm(r_ab + periodic_vectors(lattice) * shift)
        real_space_energy += erfc(sqrt(alpha) * r) / (2 * r)
    end
    return real_space_energy
end

function reciprocal_space_sum(
        depth::AbstractVector{Int},
        frac_pos_a::AbstractVector{T},
        frac_pos_b::AbstractVector{T},
        lattice::Lattice{T},
        alpha::T
        ) where T
    
end