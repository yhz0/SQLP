"""
A representation of a cut under an epigrpah in SD.
Stands for η≥α+βx.
"""
struct sdCut
    # Coefficients stands for alpha + beta' x
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
    # Epigraph weights relative to the master problem.
    epigraph_weight::Float64

    # Subproblem template coefficients shared by all scenarios.
    subproblem_coef::sdSubprobCoefficients

    # Scenario weights
    scenario_weight::Vector{Float64}
    # Total scenario weights
    total_scenario_weight::Float64

    # Scenarios
    scenario_list::Vector{spSmpsScenario}
    scenario_delta::Vector{sdDeltaCoefficients}

    sdEpigraph() = new()
end
