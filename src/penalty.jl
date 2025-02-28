function minimum_distance(ion_a::Ion{T}, ion_b::Ion{T}, lattice::Lattice{T}) where T<:Real
    return min([norm(lattice.vectors * (ion_a.frac_pos - ion_b.frac_pos + shift)) for shift in build_shifts([1, 1, 1])]...)
end

function radii_penalty(ion_a::Ion{T}, ion_b::Ion{T}, lattice::Lattice{T}, c::Float64, penalty::Float64=1e3) where T<:Real
    if ion_a ≈ ion_b
        return 0.0
    else
        return minimum_distance(ion_a, ion_b, lattice) / (ion_a.radii + ion_b.radii) > c ? 0.0 : penalty
    end
end