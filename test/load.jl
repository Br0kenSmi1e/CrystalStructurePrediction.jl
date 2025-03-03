using JuMP, Gurobi

model = read_from_file("/Users/longli/pycode/ipcsp/ipcsp/model.mps")
set_optimizer(model, Gurobi.Optimizer; add_bridges=false)
optimize!(model)
assert_is_solved_and_feasible(model)