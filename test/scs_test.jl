include("../src/TwoSD.jl")
# Driver file for scs algorithm
using .TwoSD, MathOptInterface, JuMP, GLPK, Test

const MOI = MathOptInterface

optimizer = GLPK.Optimizer

cor = TwoSD.read_cor(joinpath("spInput", "lands", "lands.cor"))
tim = TwoSD.read_tim(joinpath("spInput", "lands", "lands.tim"))
sto = TwoSD.read_sto(joinpath("spInput", "lands", "lands.sto"))

# the stage templates
sp1 = TwoSD.get_smps_stage_template(cor, tim, 1)
sp2 = TwoSD.get_smps_stage_template(cor, tim, 2)
set_optimizer(sp2.model, optimizer)

# generate template coefficients
sp1_coef = TwoSD.extract_coefficients(sp1)
sp2_coef = TwoSD.extract_coefficients(sp2)

# Initial feasible solution
x = [2., 3, 2, 5]
TwoSD.check_first_stage_feasible(sp1, x, optimizer=optimizer)

# Set up sample set (delta set) and dual vertex set
# Delta means the difference between template and scenario LPs
samples = TwoSD.spSmpsScenario[]
deltas = TwoSD.sdDeltaCoefficients[]
dv = TwoSD.sdDualVertexSet()

# generate a new scenario and add it to the scenario set
new_sample = rand(sto)
new_delta = TwoSD.delta_coefficients(sp2_coef, new_sample)
push!(samples, new_sample)
push!(deltas, new_delta)

# solve subproblem corresponding to that scenario
# record dual
obj, y_opt, dual_opt = TwoSD.solve_problem!(sp2, x, new_sample)
push!(dv, dual_opt)

# argmax procedure computes for each scenario the max value and the dual that takes the value
max_vals, max_args = TwoSD.argmax_procedure(sp2_coef, deltas, x, dv)
