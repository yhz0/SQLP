using SQLP, MathOptInterface, JuMP, CPLEX, Test

optimizer = CPLEX.Optimizer
cor = SQLP.read_cor(joinpath("spInput", "lands", "lands.cor"))
tim = SQLP.read_tim(joinpath("spInput", "lands", "lands.tim"))
sto = SQLP.read_sto(joinpath("spInput", "lands", "lands.sto"))

sp1 = SQLP.get_smps_stage_template(cor, tim, 1)
sp2 = SQLP.get_smps_stage_template(cor, tim, 2)

cell = SQLP.sdCell(sp1)
epi1 = SQLP.sdEpigraph(sp2, 0.5, 0.0)
epi2 = SQLP.sdEpigraph(sp2, 0.5, 0.0)
SQLP.bind_epigraph!(cell, epi1)
SQLP.bind_epigraph!(cell, epi2)

set_optimizer(cell.master, optimizer)
set_optimizer(epi1.prob.model, optimizer)
set_optimizer(epi2.prob.model, optimizer)

# A starting solution
x0 = [3., 3, 3, 3]
@test SQLP.check_first_stage_feasible(sp1, x0; optimizer)
cell.x_incumbent .= x0
cell.x_candidate .= x0

# TODO debug this procedure
SQLP.standard_sd_iteration!(cell, [rand(sto), rand(sto)])

epi1.cuts

cell

