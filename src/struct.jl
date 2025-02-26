struct Lattice{T}
    vectors::Matrix{T}
    pbc::Vector{Bool}
end

struct Ion{T}
    species::Symbol
    charge::Int
    radii::Float64
    frac_pos::Vector{T}
end

periodic_vectors(lattice::Lattice) = lattice.vectors[:, lattice.pbc]
build_grid(nsize::AbstractVector{Int}) = map(x->collect(x.I .- 1), vec(CartesianIndices(Tuple(nsize))))