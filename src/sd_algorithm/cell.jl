"""
Master problem and intermediate solver information.
"""
mutable struct sdCell
    # Master Problem
    master::Model

    # Reference to first stage variables, constraints, epigrpah in master problem
    x_ref::Vector{VariableRef}
    root_stage_con::Vector{ConstraintRef}
    epivar_ref::Vector{VariableRef}

    # Epigraph information
    epi::Vector{sdEpigraph}

    # Regularization strength
    reg::Float64

    sdCell() = new()
end

"""
Copy master problem and master constraints into the cell.
"""
function initialize_cell!(cell::sdCell, sp1::spStageProblem)
    cell.master, refmap = JuMP.copy_model(sp1.model)
    cell.x_ref = [refmap[v] for v in sp1.current_stage_vars]
    cell.root_stage_con = [refmap[c] for c in sp1.stage_constraints]
    cell.epivar_ref = []
    return
end

"""
Add a scenario to epigraph variable epi.
"""
function add_scenario(epi::sdEpigraph, scenario::spSmpsScenario, weight::Float64)
    push!(epi.scenarios, scenario)
    push!(epi.scenario_weight, weight)
    epi.total_scenario_weight += weight
end
