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

"""
A representation of a cut under an epigrpah in SD.
Stands for η≥α+βx.
"""
struct sdCut
    alpha::Float64
    beta::Vector{Float64}
    # Records when this cut is generated. Used to scale the cut as appropriate.
    weight_mark::Float64
    incumbent::Bool
end

"""
A convex piecewise function approximation.
"""
mutable struct sdEpigraph
    # Epigraph weights
    epigraph_weight::Float64

    # Scenario weights
    scenario_weight::Vector{Float64}

    # Total scenario weights
    total_scenario_weight::Float64

    # Scenarios
    scenarios::Vector{Any}

    sdEpigraph() = new()
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

"""
Add a scenario to epigraph variable epi.
"""
function add_scenario(epi::sdEpigraph, scenario, weight::Float64)
    push!(epi.scenarios, scenario)
    push!(epi.scenario_weight, weight)
    epi.total_scenario_weight += weight
end
