# Implementation of Two-stage Stochastic Decomposition algorithm

using LinearAlgebra

include("subprob.jl")
include("epigraph.jl")
include("cell.jl")

"""
Do an iteration of SD given the scenario list.
"""
function standard_sd_iteration(cell::sdCell, scenario_list::Vector{spSmpsScenario})
    # Make sure the length of the epigraph variable is the same.
    @assert(length(scenario_list) == length(cell.epivar_ref))

    error("Unimplemented")
end