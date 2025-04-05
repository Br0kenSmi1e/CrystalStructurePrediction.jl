module CrystalStructurePrediction

using SpecialFunctions
using LinearAlgebra
using JuMP
using SCIP
using StaticArrays

export Lattice, Ion, IonType
export ions_on_grid, interaction_energy
export optimize_linear, optimize_quadratic

include("struct.jl")
include("interaction.jl")
include("build_problem.jl")

end
