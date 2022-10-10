using JuMP

"""
The problem contains oracle for which SD operates on.
It should support the following
"""
abstract type spProblem end

"""
`solve_subproblem(prob, in_state, obs)` solves the subproblem given input state
and observation. Function should return a structure (TODO) giving the solution
of the subproblem.
"""
function solve_subproblem(prob::spProblem, in_state, obs)
    error("unimplemented")
end

"""
`evaluate_subproblem(prob, in_state, dual_point, obs)` evaluate the subproblem
at a given input state and dual point, and an observation.
"""
function evaluate_subproblem(prob::spProblem, instate)
    error("unimplemented")
end
