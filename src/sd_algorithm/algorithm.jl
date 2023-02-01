"""
Do an iteration of SD given the scenario list.
One SD iteration proceeds as follows:
0. A list of new scenario is provided for each epigraph.
Initialized cell is also provided.
1. Add all the scenarios to their assigned epigraph variable.
Do 2, 3 (in parallel) for each epigraph variable.
2. Solve the newly added scenario subproblems within the assigned
epigraph variable, at candidate solution. Add the dual vertices
to cell's dual vertex Set. (Need mutex on DV set)
3. Do 2, with candidate replaced by incumbent. (Need mutex on DV set)
4. If the cell's master problem was solved before, go to each epigraph
and remove cut that has positive dual multipliers.
5. In parallel, generate new candidate cuts from each epigraph variable
by argmax procedure, at the candidate solution. Store it in the epigraph.
6. Do 5, with incumbent solution.
7. Sync the cuts to the cell's master. Book-keeping the constraints.
Solve master to get new candidate solution.
8. Do incumbent selection. Replace incumbent if needed.
"""
function standard_sd_iteration!(cell::sdCell, scenario_list::Vector{spSmpsScenario})
    # Make sure the length of the epigraph variable is the same.
    @assert(length(scenario_list) == length(cell.epivar_ref))

    for i in eachindex(scenario_list)
        add_scenario!(cell.epi[i], scenario_list[i], 1.0)
    end

    # TODO
    error("Unimplemented")
end
