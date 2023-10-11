using .TwoSD, JuMP, GLPK, Test

optimizer = GLPK.Optimizer
cor = TwoSD.read_cor(joinpath("spInput", "lands", "lands.cor"))
tim = TwoSD.read_tim(joinpath("spInput", "lands", "lands.tim"))
sto = TwoSD.read_sto(joinpath("spInput", "lands", "lands.sto"))

# the stage templates
sp1 = TwoSD.get_smps_stage_template(cor, tim, 1)
sp2 = TwoSD.get_smps_stage_template(cor, tim, 2)

# show model templates
# print(sp1.model)
# print(sp2.model)
set_optimizer(sp2.model, optimizer)

# Generate a scenario and instantiate
sc = [TwoSD.spSmpsPosition("RHS", "S2C5") => 7.0]
TwoSD.instantiate!(sp2, sc)

# solve subproblem and get subgradient
X = [2.0, 3.0, 4.0, 5.0]
fix.(sp2.last_stage_vars, X; force=true)

optimize!(sp2.model)
subgrad = dual.(FixRef.(sp2.last_stage_vars))
@test length(subgrad) == 4
@test subgrad == [-11.0, -6.0, -19.0, 0.0]
