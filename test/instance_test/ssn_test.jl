using .TwoSD, MathOptInterface, JuMP, CPLEX, Test, Printf

Base.show(io::IO, f::Float64) = @printf(io, "%.4f", f)

optimizer = CPLEX.Optimizer
# Load files
prob_name = "ssn"
cor = TwoSD.read_cor(joinpath("spInput", "$prob_name", "$prob_name.cor"))
tim = TwoSD.read_tim(joinpath("spInput", "$prob_name", "$prob_name.tim"))
sto = TwoSD.read_sto(joinpath("spInput", "$prob_name", "$prob_name.sto"))

sp1 = TwoSD.get_smps_stage_template(cor, tim, 1)
sp2 = TwoSD.get_smps_stage_template(cor, tim, 2)

set_optimizer(sp2.model, optimizer)

# Set up cell
cell = TwoSD.sdCell(sp1)
epi1 = TwoSD.sdEpigraph(sp2, 1.0, 0.0)
TwoSD.bind_epigraph!(cell, epi1)

set_optimizer(cell.master, optimizer)
set_silent(cell.master)
for epi in cell.epi
    set_optimizer(epi.prob.model, optimizer)
    set_silent(epi.prob.model)
end

# A starting solution
# ssn
x0 = zeros(89)
@test TwoSD.check_first_stage_feasible(sp1, x0; optimizer)
cell.x_incumbent .= x0
cell.x_candidate .= x0

# Populate with initial samples
# for i = 1:1000
#     TwoSD.add_scenario!(cell.epi[1], rand(sto))
# end

using Random

function test()
    Random.seed!(42)
    qss = TwoSD.AdaptiveQuadScalarSchedule()
    cell.ext[:quad_scalar]=0.001
    for i = 1:3000
        TwoSD.sd_iteration!(cell, [rand(sto)]; quad_scalar_schedule=qss)

        lb = cell.improvement_info.candidate_estimation
        repl = cell.improvement_info.is_improved
        # ub = NaN
        # if i % 100 == 0
        #     ub = TwoSD.evaluate(sp1, sp2, sto, cell.x_incumbent; N=10000)
        # else
        #     ub = NaN
        # end
        
        println("it $i quad $(cell.ext[:quad_scalar])")
    end
end
test()
# @enter test()