using CPLEX, JuMP

model = direct_model(optimizer_with_attributes(
    CPLEX.Optimizer, CPLEX.PassNames() => true
))

@variables(model, begin
    0<=x<=10
    y1 >= 0
    y2 >= 0
    y3 >= 0
end)

@objective(model, Min, 0.5*y1 + 4.0*y2)
@constraints(model, begin
    con1, y1 >= x - y3
    con2, y2 >= y3 - x
    con3, y3 == 5
end)

bm = backend(model)
CPLEX.CPXwriteprob(bm.env, bm.lp, "newsvendor.mps", "MPS")