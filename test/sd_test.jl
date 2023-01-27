using SQLP, JuMP, CPLEX, Test
cor = SQLP.read_cor(joinpath("spInput", "lands", "lands.cor"))
tim = SQLP.read_tim(joinpath("spInput", "lands", "lands.tim"))
sto = SQLP.read_sto(joinpath("spInput", "lands", "lands.sto"))

# the stage templates
sp1 = SQLP.get_smps_stage_template(cor, tim, 1)
sp2 = SQLP.get_smps_stage_template(cor, tim, 2)
set_optimizer(sp2.model, CPLEX.Optimizer)
set_silent(sp2.model)

# Create cell
cell = SQLP.sdCell()
@test cell !== nothing

# Initialze master
SQLP.initialize_master!(cell, sp1)

@test length(cell.x_ref) == 4
@test length(cell.root_stage_con) == 2
@test cell.x_ref[1].model === cell.master

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
