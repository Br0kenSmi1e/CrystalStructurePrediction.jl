struct Lattice{D, T, L}
    vectors::SMatrix{D, D, T, L}
    pbc::NTuple{D, Bool}
end
function Lattice(vectors::AbstractMatrix{T}, pbc::NTuple{D, Bool}) where {D, T}
    return Lattice(SMatrix{D, D}(vectors), pbc)
end

struct Ion{D, T}
    species::Symbol
    charge::Int
    radii::Float64
    frac_pos::SVector{D, T}
end
function Ion(species::Symbol, charge::Int, radii::Float64, frac_pos)
    return Ion(species, charge, radii, SVector(frac_pos...))
end

periodic_vectors(lattice::Lattice) = lattice.vectors[:, findall(lattice.pbc)]
build_grid(nsize::NTuple{N, Int}) where N = map(x->SVector{N, Int}(x.I .- 1), vec(CartesianIndices(nsize)))

import Base.isapprox
Base.isapprox(a::Ion, b::Ion) = (a.species == b.species) && (a.charge == b.charge) && (a.radii ≈ b.radii) && (norm(a.frac_pos - b.frac_pos - floor.(a.frac_pos - b.frac_pos))≈0)