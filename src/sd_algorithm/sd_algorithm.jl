using LinearAlgebra

include("subprob.jl")

# Implementation of Two-stage Stochastic Decomposition algorithm

mutable struct sdCell
    # Master Problem
    master::Model

    # Reference to first stage variables
    x_ref::Vector{VariableRef}

    # Reference to first stage constraints
    root_stage_con::Vector{ConstraintRef}

    # Reference to epigraph variables
    epivar_ref::Vector{VariableRef}

    # Regularization strength
    reg::Float64

    sdCell() = new()
end

mutable struct sdEpigraph
    # Epigraph weights
    epigraph_weight::Float64

    # Scenario weights
    scenario_weight::Vector{Float64}

    # Scenarios
    scenarios::Vector{Any}
end


"""
Copy master problem and master constraints into the cell.
"""
function initialize_master!(cell::sdCell, root_prob::spStageProblem)
    cell.master, refmap = JuMP.copy_model(root_prob.model)

    cell.x_ref = [refmap[v] for v in root_prob.current_stage_vars]
    cell.root_stage_con = [refmap[c] for c in root_prob.stage_constraints]

    return
end

