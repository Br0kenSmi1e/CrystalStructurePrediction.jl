using CrystalStructurePrediction

grid_size = [4, 4, 4]
population_list = [1, 1, 3]
species_list = ["Sr", "Ti", "O"]
charge_list = [+2, +4, -2]
L = [3.94513 0 0; 0 3.94513 0; 0 0 3.94513]
depth = [4, 4, 4]
alpha = 2.0 / 3.94513

lattice = Lattice(L, [true, true, true])
ion_list = build_ion_list(grid_size, species_list, charge_list)
matrix = build_matrix(ion_list, lattice, alpha, depth, depth, depth)
for i in range(1,prod(grid_size)*length(population_list))
for j in range(1,prod(grid_size)*length(population_list))
    if isnan(matrix[i,j])
        @show i,j
        break
        break
    end
end
end
# @show matrix[33, 45]
solution = build_problem(grid_size, population_list, matrix)
@show solution