using SQLP, MathOptInterface, JuMP, Test, Printf


Base.show(io::IO, f::Float64) = @printf(io, "%.4f", f)

using CPLEX, Ipopt
master_optimizer = CPLEX.Optimizer
subprob_optimizer = CPLEX.Optimizer

# Load files
prob_name = "baa99-20"
cor = SQLP.read_cor(joinpath("spInput", "$prob_name", "$prob_name.cor"))
tim = SQLP.read_tim(joinpath("spInput", "$prob_name", "$prob_name.tim"))
sto = SQLP.read_sto(joinpath("spInput", "$prob_name", "$prob_name.sto"))

sp1 = SQLP.get_smps_stage_template(cor, tim, 1)
sp2 = SQLP.get_smps_stage_template(cor, tim, 2)

set_optimizer(sp2.model, subprob_optimizer)

# Set up cell
cell = SQLP.sdCell(sp1)
epi1 = SQLP.sdEpigraph(sp2, 1.0, -500000.0)
SQLP.bind_epigraph!(cell, epi1)

set_optimizer(cell.master, master_optimizer)
set_silent(cell.master)
for epi in cell.epi
    set_optimizer(epi.prob.model, subprob_optimizer)
    set_silent(epi.prob.model)
end

# A starting solution
# lands
# x0 = [3., 3, 3, 3]
# transship
# x0 = [100.0, 200.0, 150.0, 170.0, 180.0, 170.0, 170.0]
# x0 = [10.0, 20.0, 15.0, 17.0, 18.0, 17.0, 17.0]
# ssn
# x0 = zeros(89)

# Choose a starting solution with 10 scenario problem
scenario_set = [rand(sto) for i = 1:10]
all_in_one_model, annotation = SQLP.all_in_one(sp1, sp2,scenario_set)
set_optimizer(all_in_one_model, CPLEX.Optimizer)
optimize!(all_in_one_model)
x0 = value.(annotation[0])

# x0 = zeros(20)

# Starting solution
start_obj = SQLP.evaluate(sp1, sp2, sto, x0; N=1000)
@show start_obj

cell.x_incumbent .= x0
cell.x_candidate .= x0

# Populate with initial samples
# for i = 1:100
#     SQLP.add_scenario!(cell.epi[1], rand(sto))
# end

using Random

function test()
    Random.seed!(42)
    # qss = SQLP.AdaptiveQuadScalarSchedule()
    # cell.ext[:quad_scalar] = 0.001
    for i = 1:1000
        SQLP.sd_iteration!(cell, [rand(sto)]; quad_scalar_schedule=SQLP.ConstantQuadScalarSchedule(0.1))

        if i % 10 == 0
            lb = cell.improvement_info.candidate_estimation
            repl = cell.improvement_info.is_improved

            if i % 100 == 0
                ub = SQLP.evaluate(sp1, sp2, sto, cell.x_incumbent; N=1000)
            else
                ub = NaN
            end

            qs = cell.ext[:quad_scalar]
            println("Iter $i lb=$lb ub=$ub repl=$repl qs=$qs dual=$(length(cell.dual_vertices))")
            @show cell
        end
    end
end

test()