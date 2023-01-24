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

# Test extract coefficients
coef = SQLP.extract_coefficients(sp2)
