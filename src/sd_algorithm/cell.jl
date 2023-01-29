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

end

"""
Copy master problem and master constraints into the cell.
"""
function sdCell(root_prob::spStageProblem)
    # Initialize master problem by copying
    new_root_prob::spStageProblem = copy(root_prob)
    master = new_root_prob.model
    x_ref = new_root_prob.current_stage_vars
    root_stage_con = new_root_prob.stage_constraints

    epivar_ref = []
    epi = []
    reg = 0.0
    return sdCell(master, x_ref, root_stage_con, epivar_ref, epi, reg)
end


"""
Do an iteration of SD given the scenario list.
"""
function standard_sd_iteration!(cell::sdCell, scenario_list::Vector{spSmpsScenario})
    # Make sure the length of the epigraph variable is the same.
    @assert(length(scenario_list) == length(cell.epivar_ref))

    for i in eachindex(scenario_list)
        add_scenario!(cell.epi[i], scenario_list[i], 1.0)
    end

    error("Unimplemented")
end

"""
Check stopping criteria.
"""
function stopping_criteria()::Bool
    error("Unimplemented")
end