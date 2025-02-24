struct Lattice{T}
    vectors::AbstractMatrix{T}
    pbc::AbstractVector{Bool}
end

struct Ion{T}
    species::String
    frac_pos::AbstractVector{T}
    charge::Int
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
