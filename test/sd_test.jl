using SQLP, JuMP, CPLEX, Test
cor = SQLP.read_cor(joinpath("spInput", "lands", "lands.cor"))
tim = SQLP.read_tim(joinpath("spInput", "lands", "lands.tim"))
sto = SQLP.read_sto(joinpath("spInput", "lands", "lands.sto"))

# the stage templates
sp1 = SQLP.get_smps_stage_template(cor, tim, 1)
sp2 = SQLP.get_smps_stage_template(cor, tim, 2)

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

# TODO: test for eval_dual
my_dual = randn(7)
my_x = randn(4)
my_scenario = [SQLP.spSmpsPosition("RHS", "S2C5") => 5.0]
my_scenario_2 = [SQLP.spSmpsPosition("RHS", "S2C5") => 3.0]
@show SQLP.eval_dual(coef, my_x, my_dual, my_scenario)
@show SQLP.eval_dual(coef, my_x, my_dual, my_scenario_2)

