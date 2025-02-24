struct Lattice{T}
    vectors::AbstractMatrix{T}
    pbc::AbstractVector{Bool}
end

struct Ion{T}
    species::String
    frac_pos::AbstractVector{T}
    charge::Int
end

function periodic_vectors(lattice::Lattice{T}) where T
    periodic_vectors = Matrix{T}(undef, size(lattice.vectors)[1], sum(lattice.pbc))
    for column in range(1, size(lattice.vectors)[2])
        if lattice.pbc[column]
            periodic_vectors[:, column] = lattice.vectors[:, column]
        end
    end
    return periodic_vectors
end

function build_grid(nsize::AbstractVector{Int})
    function decompose(n::Int)
        remain = n
        result = Vector{Int}(undef, length(nsize))
        for i in range(1, length(nsize))
            result[i] = remain ÷ prod(nsize[i+1:end])
            remain = remain % prod(nsize[i+1:end])
        end
        return result
    end
    return [decompose(n) for n in range(0, prod(nsize) - 1)]
end