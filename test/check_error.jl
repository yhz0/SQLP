using JuMP, CPLEX
model = read_from_file("error_model.mof.json")
set_optimizer(model, CPLEX.Optimizer)
print(model)
optimize!(model)
