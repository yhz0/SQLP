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
function stopping_criteria(cell::sdCell)::Bool
    error("Unimplemented")
end