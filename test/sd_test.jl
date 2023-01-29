using SQLP, MathOptInterface, JuMP, GLPK, Test
const MOI = MathOptInterface

optimizer = GLPK.Optimizer

cor = SQLP.read_cor(joinpath("spInput", "lands", "lands.cor"))
tim = SQLP.read_tim(joinpath("spInput", "lands", "lands.tim"))
sto = SQLP.read_sto(joinpath("spInput", "lands", "lands.sto"))

# the stage templates
sp1 = SQLP.get_smps_stage_template(cor, tim, 1)
sp2 = SQLP.get_smps_stage_template(cor, tim, 2)
set_optimizer(sp2.model, optimizer)

# ===== subprob.jl
# Test extract coefficients
coef = SQLP.extract_coefficients(sp2)

@test coef.row_lookup["S2C5"] == 5
@test coef.col_lookup["X2"] == 2

# Make sure only the existing key can be found
@test_throws KeyError coef.col_lookup["Y11"]

# Test modify coefficients
my_scenario = [SQLP.spSmpsPosition("RHS", "S2C5") => 5.0]
my_scenario_2 = [SQLP.spSmpsPosition("RHS", "S2C5") => 3.0]
my_scenario_3 = [SQLP.spSmpsPosition("RHS", "S2C5") => 7.0]

SQLP.modify_coefficients!(coef, my_scenario)
@test coef.rhs[coef.row_lookup["S2C5"]] == 5.0
SQLP.modify_coefficients!(coef, my_scenario_2)
@test coef.rhs[coef.row_lookup["S2C5"]] == 3.0

# Test for delta coefficients
SQLP.instantiate!(sp2, my_scenario_2) # RHS, S2C5 => 3.0
coef = SQLP.extract_coefficients(sp2)
delta = SQLP.delta_coefficients(coef, my_scenario) # RHS, S2C5 => 5.0, so delta for RHS, S2C5 == 2.0
ind = coef.row_lookup["S2C5"]
@test delta.delta_rhs[ind] == 2.0
@test sum(delta.delta_transfer) == 0.0

# Test for eval_dual
# Calculate true answers first.
my_x = [3.0, 3.0, 3.0, 3.0]

fix.(sp2.last_stage_vars, my_x, force=true)
SQLP.instantiate!(sp2, my_scenario)
optimize!(sp2.model)
my_dual = dual.(sp2.stage_constraints)
sc_val1_true = objective_value(sp2.model)

SQLP.instantiate!(sp2, my_scenario_2)
optimize!(sp2.model)
my_dual_2 = dual.(sp2.stage_constraints)
sc_val2_true = objective_value(sp2.model)

# Now for eval_dual
delta = SQLP.delta_coefficients(coef, my_scenario)
delta2 = SQLP.delta_coefficients(coef, my_scenario_2)

sc_val1 = SQLP.eval_dual(coef, delta, my_x, my_dual)
sc_val2 = SQLP.eval_dual(coef, delta2, my_x, my_dual_2)
@test sc_val1_true == sc_val1
@test sc_val2_true == sc_val2

# Test argmax procedure
# Two solutions
x1 = [3.0, 3.0, 3.0, 3.0]
x2 = [2.0, 4.0, 2.0, 6.0]
@test SQLP.check_first_stage_feasible(sp1, x1; optimizer=optimizer)
@test SQLP.check_first_stage_feasible(sp1, x2; optimizer=optimizer)

# Some scenarios
coef = SQLP.extract_coefficients(sp2)
scenario_set = SQLP.spSmpsScenario[my_scenario, my_scenario, my_scenario_2, my_scenario_3]
delta_set = SQLP.delta_coefficients.(Ref(coef), scenario_set)

# Generate some dual points with x1
dual_points = Set{Vector{Float64}}()
for scenario in scenario_set
    v, y, p = SQLP.solve_problem!(sp2, x1, scenario)
    push!(dual_points, p)
end
@test length(dual_points) == 3

val, arg = SQLP.argmax_procedure(coef, delta_set, x2, dual_points)

# These dual points are sufficient.
# so the objective value should be the same as solving it directly.
for i in eachindex(scenario_set)
    v, y, p = SQLP.solve_problem!(sp2, x2, scenario_set[i])
    @test val[i] == v
end

# Another test case
begin
    x = [2, 3, 4, 5.]
    v, y, p = SQLP.solve_problem!(sp2, x, my_scenario_3)
    d = SQLP.delta_coefficients(coef, my_scenario_3)
    subgrad = - (d.delta_transfer + coef.transfer)' * p
    @test subgrad == [-11.0, -6.0, -19.0, 0.0]
end

# === cell.jl
cell = SQLP.sdCell(sp1)
@test cell.master !== nothing
# Test if the references are copied correctly
@test cell.master !== sp1.model
@test is_valid(cell.master, cell.x_ref[1])
@test !is_valid(cell.master, sp1.current_stage_vars[1])


# === epigraph.jl
epi = SQLP.sdEpigraph(sp2, 0.5)
# Test copy
@test epi.prob !== nothing
@test epi.prob.model !== sp2.model
@test is_valid(epi.prob.model, epi.prob.current_stage_vars[1])
@test !is_valid(epi.prob.model, sp2.current_stage_vars[1])

epi2 = SQLP.sdEpigraph(sp2, 0.5)
push!(cell.epi, epi)
push!(cell.epi, epi2)
