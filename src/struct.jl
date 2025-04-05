"""
    Lattice{D, T, L}
    Lattice(vectors::AbstractMatrix{T}, pbc::NTuple{D, Bool}) where {D, T}

A lattice is a set of vectors that define the unit cell of a crystal.

# Fields
- `vectors::SMatrix{D, D, T, L}`: The vectors that define the unit cell of a crystal.
- `pbc::NTuple{D, Bool}`: The periodic boundary conditions.
"""
struct Lattice{D, T, L}
    vectors::SMatrix{D, D, T, L}
    pbc::NTuple{D, Bool}
end
function Lattice(vectors::AbstractMatrix{T}, pbc::NTuple{D, Bool}) where {D, T}
    return Lattice(SMatrix{D, D}(vectors), pbc)
end

"""
    Ion{D, T}
    Ion(species::Symbol, charge::Int, radii::Float64, frac_pos)

An ion is a species with a charge, a radius, and a fractional position.
`D` is the dimension of the space, and `T` is the type of the coordinates.

# Fields
- `species::Symbol`: The species of the ion.
- `charge::Int`: The charge of the ion.
- `radii::Float64`: The radius of the ion.
- `frac_pos::SVector{D, T}`: The fractional position of the ion.

# Examples
To define an ion with species "O", charge -2, radius 1.35, and fractional position [0.25, 0.25, 0.25], run
```julia
julia> Ion(:O, -2, 1.35, [0.25, 0.25, 0.25])
Ion{3, Float64}(:O, -2, 1.35, [0.25, 0.25, 0.25])
```
"""
struct Ion{D, T}
    species::Symbol
    charge::Int
    radii::Float64
    frac_pos::SVector{D, T}
end
function Ion(species::Symbol, charge::Int, radii::Float64, frac_pos)
    return Ion(species, charge, radii, SVector(frac_pos...))
end

# filter out the periodic unit vectors
periodic_vectors(lattice::Lattice) = lattice.vectors[:, findall(lattice.pbc)]

Base.isapprox(a::Ion, b::Ion) = (a.species == b.species) && (a.charge == b.charge) && (a.radii ≈ b.radii) && (norm(a.frac_pos - b.frac_pos - floor.(a.frac_pos - b.frac_pos))≈0)