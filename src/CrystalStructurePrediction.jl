module CrystalStructurePrediction

using SpecialFunctions
using LinearAlgebra
using JuMP
using HiGHS
using StaticArrays

export Lattice, Ion
export build_ion_list, build_vector
export build_matrix, interaction_energy
export build_linear_problem, build_quadratic_problem

include("struct.jl")
include("interaction.jl")
include("build_matrix.jl")
include("build_problem.jl")

end
