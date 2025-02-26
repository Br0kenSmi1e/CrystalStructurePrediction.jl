using CrystalStructurePrediction

grid_size = [2, 2, 2]
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
alpha = 1 / 3.899

lattice = Lattice(L, [true, true, true])
ion_sheet = build_ion_sheet(grid_size, species_list, charge_list, radii_list)
matrix = build_matrix(ion_sheet, lattice, interaction_energy, (alpha, depth, depth, depth))
# matrix = build_matrix(ion_sheet, lattice, buckingham_sum, Tuple([depth]))
# @show matrix

energy, solution_x, solution_s = build_linear_problem(grid_size, population_list, matrix)
for index in CartesianIndices(ion_sheet)
    if solution_x[index] ≈ 1
        @show ion_sheet[index]
    end
end

# @show solution_s

@show energy

# solution = build_quadratic_problem(grid_size, population_list, matrix)
