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
    # Q: why pbc is not used?
    pbc::NTuple{D, Bool}
end
function Lattice(vectors::AbstractMatrix{T}, pbc::NTuple{D, Bool}) where {D, T}
    return Lattice(SMatrix{D, D}(vectors), pbc)
end
# convert fractional coordinates to Cartesian coordinates
cartesian(lt::Lattice, v) = lt.vectors * v
# convert Cartesian coordinates to fractional coordinates
fractional(lt::Lattice, v) = lt.vectors \ v
# reciprocal lattice vectors
reciprocal_vectors(lt::Lattice) = 2π .* transpose(inv(lt.vectors))
# volume of the unit cell
volume(lt::Lattice) = abs(det(lt.vectors))

# the minimum distance between two ions in a lattice
# TODO: this is the notoriously hard closest vector problem, try solve it with maybe integer programming?
function minimum_distance(frac_pos_a::AbstractVector{T}, frac_pos_b::AbstractVector{T}, lattice::Lattice{D, T}) where {D, T}
    return minimum(shift -> distance(lattice, frac_pos_b + SVector(shift), frac_pos_a), Iterators.product(ntuple(x->-1:1, D)...))
end
distance(a::AbstractVector{T}, b::AbstractVector{T}) where T = norm(a - b)
distance(lt::Lattice{D, T}, frac_pos_a::AbstractVector{T}, frac_pos_b::AbstractVector{T}) where {D, T} = norm(cartesian(lt, frac_pos_b - frac_pos_a))

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

"""
    IonType{T}
    IonType(species::Symbol, charge::Int, radii::T)

An ion type is a species with a charge and a radius.
`T` is the type of the coordinates.

# Fields
- `species::Symbol`: The species of the ion.
- `charge::Int`: The charge of the ion.
- `radii::T`: The radius of the ion.

# Examples
To define an ion type with species "O", charge -2, radius 1.35, run
```jldoctest
julia> IonType(:O, -2, 1.35)
IonType{Float64}(:O, -2, 1.35)
```

"""
struct IonType{T}
    species::Symbol
    charge::Int
    radii::T
end
charge(ion::IonType) = ion.charge
radii(ion::IonType) = ion.radii
species(ion::IonType) = ion.species

"""
    Ion{D, T}
    Ion(type::IonType{T}, frac_pos)

An ion is a type with a fractional position.
`D` is the dimension of the space, and `T` is the type of the coordinates.

# Fields
- `type::IonType{T}`: The type of the ion.
- `frac_pos::SVector{D, T}`: The fractional position of the ion.
"""
struct Ion{D, T}
    type::IonType{T}
    frac_pos::SVector{D, T}
end
function Ion(type::IonType{T}, frac_pos) where T
    return Ion(type, SVector(frac_pos...))
end

charge(ion::Ion) = ion.type.charge
radii(ion::Ion) = ion.type.radii
species(ion::Ion) = ion.type.species
Base.isapprox(a::Ion, b::Ion; kwargs...) = (a.type == b.type) && isapprox(a.frac_pos, b.frac_pos; kwargs...)