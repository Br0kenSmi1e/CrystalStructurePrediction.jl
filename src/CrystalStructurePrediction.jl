module CrystalStructurePrediction

using SpecialFunctions
using LinearAlgebra
using JuMP
using SCIP
using StaticArrays

export Lattice, Ion, IonType
export ions_on_grid, build_vector, build_matrix, interaction_energy
export build_linear_problem, build_quadratic_problem, build_proximal_pairs

include("struct.jl")
include("interaction.jl")
include("build_problem.jl")

end
