using CrystalStructurePrediction

grid_size = [4, 4, 4]
population_list = [1, 1, 3]
species_list = [:Sr, :Ti, :O]
charge_list = [+2, +4, -2]
radii_list = [1.18, 0.42, 1.35]
# population_list = [1, 1]
# species_list = [:Na, :Cl]
# charge_list = [+1, -1]
# radii_list = [1.18, 0.42]
L = 3.899 * [1 0 0; 0 1 0; 0 0 1]
depth = [4, 4, 4]
alpha = 2 / 3.899

lattice = Lattice(L, [true, true, true])
ion_list = build_ion_list(grid_size, species_list, charge_list, radii_list)
proximal_pairs = build_proximal_pairs(ion_list, lattice, 0.75)

# vector = build_vector(ion_list, lattice, interaction_energy, (alpha, depth, depth, depth))
# energy, solution_x, solution_s = build_linear_problem(grid_size, population_list, vector, proximal_pairs)

matrix = build_matrix(ion_list, lattice, interaction_energy, (alpha, depth, depth, depth))
# csp = build_quadratic_problem(grid_size, population_list, matrix, proximal_pairs)
# @show csp
# @show matrix
energy, solution_x, csp = build_quadratic_problem(grid_size, population_list, matrix, proximal_pairs)
for index in CartesianIndices(ion_list)
    if solution_x[index] ≈ 1
        @show ion_list[index]
    end
end

# # @show solution_s

# @show energy

# solution = build_quadratic_problem(grid_size, population_list, matrix)
