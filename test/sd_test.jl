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
    local v, y, p = SQLP.solve_problem!(sp2, x1, scenario)
    push!(dual_points, p)
end
@test length(dual_points) == 3

val, arg = SQLP.argmax_procedure(coef, delta_set, x2, dual_points)

# These dual points are sufficient.
# so the objective value should be the same as solving it directly.
for i in eachindex(scenario_set)
    local v, y, p = SQLP.solve_problem!(sp2, x2, scenario_set[i])
    @test val[i] == v
end

# Another test case
begin
    x = [2, 3, 4, 5.]
    local v, y, p = SQLP.solve_problem!(sp2, x, my_scenario_3)
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

@test length(cell.x_candidate) == 4
@test cell.objf == QuadExpr([10, 7, 16, 6]' * cell.x_ref)

# === epigraph.jl
epi = SQLP.sdEpigraph(sp2, 0.5, 0.0)
# Test copy
@test epi.prob !== nothing
@test epi.prob.model !== sp2.model
@test is_valid(epi.prob.model, epi.prob.current_stage_vars[1])
@test !is_valid(epi.prob.model, sp2.current_stage_vars[1])

# Epigraph With lower bounds
epi2 = SQLP.sdEpigraph(sp2, 0.5, 100.0)

# Test binding epigraph
SQLP.bind_epigraph!(cell, epi)
SQLP.bind_epigraph!(cell, epi2)
@test epi === cell.epi[1]
@test length(cell.epi) == 2
@test length(cell.epicon_ref) == 2
@test length(cell.epicon_incumbent_ref) == 2
@test length(cell.epivar_ref) == 2
@test is_valid(cell.master, cell.epivar_ref[1])
@test lower_bound(cell.epivar_ref[1]) == cell.epi[1].lower_bound

@test cell.objf == 
    QuadExpr([10, 7, 16, 6]' * cell.x_ref + [0.5, 0.5]' * cell.epivar_ref)

# Add scenario to epigraph
SQLP.add_scenario!(cell.epi[1], my_scenario, 1.0)
SQLP.add_scenario!(cell.epi[1], my_scenario_2, 1.0)

SQLP.add_scenario!(cell.epi[2], my_scenario_2, 1.0)
SQLP.add_scenario!(cell.epi[2], my_scenario_3, 1.0)

@test cell.epi[1].total_scenario_weight == 2.0

# === Cuts related
# Test add_cut_to_master!
# cut eta >= 1 + [2, 3, 4, 5] x with created at weight = 0.1
cut = SQLP.sdCut(1.0, [2.0, 3.0, 4.0, 5.0], 0.1)
con = SQLP.add_cut_to_master!(cell.master, cut, cell.epivar_ref[1], cell.x_ref,
    1.0, 0.0)
@test constraint_object(con).func == cell.epivar_ref[1] - [2, 3, 4, 5]' * cell.x_ref

# Test remove_cuts!
# Register the constraint as test case then remove it.
push!(cell.epicon_ref[1], con)
SQLP.remove_cuts!(cell, 1)
@test !is_valid(cell.master, con)
@test isempty(cell.epicon_ref[1])

# Test sync_cuts!
# Two demo cuts, and an incumbent cut
cut1 = SQLP.sdCut(1.0, [2.0, 3.0, 4.0, 5.0], 1.0)
cut2 = SQLP.sdCut(6.0, [7.0, 8.0, 9.0, 10.0], 2.0)
inc_cut = SQLP.sdCut(11, [12.0, 13.0, 14.0, 15.0], 1.0)
push!(cell.epi[1].cuts, cut1)
push!(cell.epi[1].cuts, cut2)
cell.epi[1].incumbent_cut = inc_cut
push!(cell.epi[2].cuts, cut1)

SQLP.sync_cuts!(cell, cell.epi[1], 1)
@test length(cell.epicon_ref[1]) == 2
@test cell.epicon_incumbent_ref[1] !== nothing

SQLP.sync_cuts!(cell)
# Test that no redundant cut is added
@test length(cell.epicon_ref[1]) == 2

# Test lower bounds on epigraphs are added correctly when discounted
# Epi2 has lb 100, and cut is discounted by 0.5. rhs was 1.0,
# so the added cut should has rhs 100*0.5+1.0*0.5=50.5
@test normalized_rhs(cell.epicon_ref[2][1]) == 50.5

# Test evaluate_epigraph
# TODO: rewrite tests
@test SQLP.evaluate_epigraph(cell.epi[1], [10.0, 10.0, 10.0, 10.0]) == 551.0*0.5
# # test if the discount is applied correctly and the lower bound is applied
@test SQLP.evaluate_epigraph(cell.epi[2], [10.0, 10.0, 10.0, 10.0]) == (141/2 + 100/2)*0.5
@test SQLP.evaluate_epigraph(cell.epi[2], [-1.0, -1, -1, -1]) == 100.0*0.5

# Test extracting sdEpigraphInfo from sdEpigraph
# Make sure that mutating/deleting sdEpigraph does not invalidate sdEpigraphInfo

epi3 = SQLP.sdEpigraph(sp2, 1.0, 1.0)
push!(epi3.cuts, cut1)
push!(epi3.cuts, cut2)
epi3.incumbent_cut = inc_cut
epi_info = SQLP.sdEpigraphInfo(epi3)
empty!(epi3.cuts)
@test !isempty(epi_info.cuts)

# Test for build_sasa_cut

x = [2., 3., 4., 5.]
dv = [my_dual, my_dual_2]
# Change some weight to see if it works correctly
epi2 = SQLP.sdEpigraph(sp2, 0.5, 100.0)
SQLP.add_scenario!(epi2, my_scenario_2, 1.5) # RHS = 3.0
SQLP.add_scenario!(epi2, my_scenario_3, 0.5) # RHS = 7.0

# Manual calculation:
# Epi2 scenario #1 (3.0): dv1=168 dv2=169 choose dv[2]
# SQLP.eval_dual(epi[2].subproblem_coef, epi[2].scenario_delta[1], x, dv[1])
# SQLP.eval_dual(epi[2].subproblem_coef, epi[2].scenario_delta[1], x, dv[2])
# Epi2 scenario #2 (7.0): dv1=344 dv2=331 choose dv[1]
# SQLP.eval_dual(epi[2].subproblem_coef, epi[2].scenario_delta[2], x, dv[1])
# SQLP.eval_dual(epi[2].subproblem_coef, epi[2].scenario_delta[2], x, dv[2])

r1 = epi2.subproblem_coef.rhs + epi2.scenario_delta[1].delta_rhs
T1 = epi2.subproblem_coef.transfer + epi2.scenario_delta[1].delta_transfer
r2 = epi2.subproblem_coef.rhs + epi2.scenario_delta[2].delta_rhs
T2 = epi2.subproblem_coef.transfer + epi2.scenario_delta[2].delta_transfer

expected_alpha = 1.5/2.0 * dv[2]' * r1 + 0.5/2.0 * dv[1]' * r2
expected_beta = 1.5/2.0 * T1' * dv[2] + 0.5/2.0 * T2' * dv[1]

cut = SQLP.build_sasa_cut(epi2, x, dv)
@test cut.alpha == expected_alpha
@test cut.beta == expected_beta
@test cut.weight_mark == 2.0
