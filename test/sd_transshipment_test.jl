using SQLP, MathOptInterface, JuMP, CPLEX, Test, Printf

Base.show(io::IO, f::Float64) = @printf(io, "%.4f", f)

optimizer = CPLEX.Optimizer
# Load files
cor = SQLP.read_cor(joinpath("spInput", "transship", "transship.cor"))
tim = SQLP.read_tim(joinpath("spInput", "transship", "transship.tim"))
sto = SQLP.read_sto(joinpath("spInput", "transship", "transship.sto"))

sp1 = SQLP.get_smps_stage_template(cor, tim, 1)
sp2 = SQLP.get_smps_stage_template(cor, tim, 2)

set_optimizer(sp2.model, optimizer)

# Set up cell
cell = SQLP.sdCell(sp1)
epi1 = SQLP.sdEpigraph(sp2, 1.0, 0.0)
SQLP.bind_epigraph!(cell, epi1)

set_optimizer(cell.master, optimizer)
set_silent(cell.master)
for epi in cell.epi
    set_optimizer(epi.prob.model, optimizer)
    set_silent(epi.prob.model)
end


# A starting solution
x0 = zeros(7)
@test SQLP.check_first_stage_feasible(sp1, x0; optimizer)
cell.x_incumbent .= x0
cell.x_candidate .= x0

# Populate with initial samples
for i = 1:1000
    SQLP.add_scenario!(cell.epi[1], rand(sto))
    # SQLP.add_scenario!(cell.epi[2], rand(sto))
end


for i = 1:1000
    x, lb, repl = SQLP.sd_iteration!(cell, [rand(sto)]; rho=0.1)
    println("Iter $i lb=$lb repl=$repl")
end
