"""
A representation of a cut under an epigrpah in SD.
Stands for η≥α+βx.
"""
struct sdCut
    # Coefficients stands for alpha + beta' x
    alpha::Float64
    beta::Vector{Float64}
    
    # Records when this cut is generated to scale the cuts.
    weight_mark::Float64
    incumbent::Bool
end

"""
A convex piecewise function approximation.
"""
mutable struct sdEpigraph
    # Subproblem corresponding to this epigraph variable.
    prob::spStageProblem

    # Epigraph weights relative to the master problem.
    objective_weight::Float64

    # Subproblem template coefficients shared by all scenarios.
    subproblem_coef::sdSubprobCoefficients

    # Scenario weights
    scenario_weight::Vector{Float64}
    # Total scenario weights
    total_scenario_weight::Float64

    # Scenarios
    scenario_list::Vector{spSmpsScenario}
    scenario_delta::Vector{sdDeltaCoefficients}

    # Store cuts
    cuts::Vector{sdCut}
    incumbent_cut::Union{Nothing, sdCut}
end

"""
Create a new epigraph with the specified weight and subproblem template.
Make sure the subprob is not reused so coefficients change does not
overwrite each other in parallel version.
"""
function sdEpigraph(
    prob::spStageProblem,
    objective_weight::Float64)
    
    new_prob = copy(prob)
    coef = extract_coefficients(new_prob)
    
    return sdEpigraph(new_prob, objective_weight,
        coef, [], 0.0, [], [], [], nothing)
end

"""
Print an epigraph variable.
"""
function Base.show(io::IO, epi::sdEpigraph)
    objw = epi.objective_weight
    sc_cnt = length(epi.scenario_list)
    sc_w = epi.total_scenario_weight
    cut_cnt = length(epi.cuts)
    if epi.incumbent_cut !== nothing
        cut_cnt = cut_cnt + 1
    end
    print(io, "sdEpigraph objw=$objw sc_cnt=$sc_cnt sc_w=$sc_w cut_cnt=$cut_cnt")
    return
end


"""
Add a scenario to epigraph variable epi.
"""
function add_scenario!(
    epi::sdEpigraph,
    scenario::spSmpsScenario,
    weight::Float64 = 1.0)

    # Record raw scenario
    push!(epi.scenario_list, scenario)
    push!(epi.scenario_weight, weight)
    epi.total_scenario_weight += weight

    # Calculate delta coefficients, store.
    delta = delta_coefficients(epi.subproblem_coef, scenario)
    push!(epi.scenario_delta, delta)
    
    return
end
