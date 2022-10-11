"""
The problem contains oracle for which SD operates on.
It should support the following
"""
abstract type spProblem end

"""
`solve_subproblem(prob, in_state, obs)` solves the subproblem given input state
and observation. Function should return a structure (#TODO) giving the solution
of the subproblem.
"""
function get_problem(prob::spProblem, stage::Int, in_state, obs) end

"""
Evaluate the stage problem of prob at a given input state in_state
and dual point, and an observation.
"""
function evaluate_dual_subproblem(prob::spProblem, stage::Int,
    in_state, dual_point, obs) end

"""
Get (sub)gradient of primal objective, ignoring the constraints.
"""
function objective_grad(prob::spProblem, stage::Int, primal_point) end

# TODO: design randomness struct
# abstract type spRandomness end

# struct smps_randomness <: spRandomness

# end