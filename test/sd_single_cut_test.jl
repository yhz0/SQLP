using SQLP, MathOptInterface, JuMP, CPLEX, Test, Printf

Base.show(io::IO, f::Float64) = @printf(io, "%.4f", f)

optimizer = CPLEX.Optimizer
# Load files
cor = SQLP.read_cor(joinpath("spInput", "lands", "lands.cor"))
tim = SQLP.read_tim(joinpath("spInput", "lands", "lands.tim"))
sto = SQLP.read_sto(joinpath("spInput", "lands", "lands.sto"))

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
x0 = [3., 3, 3, 3]
@test SQLP.check_first_stage_feasible(sp1, x0; optimizer)
cell.x_incumbent .= x0
cell.x_candidate .= x0

# Populate with initial samples
# for i = 1:1000
#     SQLP.add_scenario!(cell.epi[1], rand(sto))
# end


hist = Float64[]
for i = 1:1000
    x, lb, repl = SQLP.sd_iteration!(cell, [rand(sto)]; rho=0.1)
    # ub = SQLP.evaluate(sp1, sp2, sto, cell.x_incumbent; N=1000)
    ub = NaN
    println("Iter $i lb=$lb ub=$ub repl=$repl")
    push!(hist, ub)
end
