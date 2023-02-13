# Implementation of main algorithm.

# One SD iteration proceeds as follows:
# 0. A list of new scenario is provided for each epigraph.
# Initialized cell is also provided.
# 1. Add all the scenarios to their assigned epigraph variable.
# Do 2, 3 (in parallel) for each epigraph variable.
# 2. Solve the newly added scenario subproblems within the assigned
# epigraph variable, at candidate solution. Add the dual vertices
# to cell's dual vertex Set. (Need mutex on DV set)
# 3. Do 2, with candidate replaced by incumbent. (Need mutex on DV set)
# 4. If the cell's master problem was solved before, go to each epigraph
# and remove cut that has positive dual multipliers.
# 5. Generate new candidate cuts from each epigraph variable
# by argmax procedure, at the candidate solution. Store it in the epigraph.
# 6. Do 5, update incumbent cut with incumbent solution.
# 7. Do incumbent selection. Replace incumbent if needed.
# 8. Sync the cuts to the cell's master. Solve master to get new candidate solution.

"""
Threshold of dual values for removing the cuts from the master program.
"""
const CUT_REMOVE_TOLERANCE = 0.001

"""
Do an iteration of SD given the scenario list.

# Arguments
- `cell`: the root cell to perform sd iteration on.
- `scenario_list::Vector{spSmpsScenario}`: a vector of new scenarios to insert to each epigraph.
Should have the same length as the number of epigraph variables in the cell.
- `quad_scalar_schedule::Function=ConstantQuadScalarSchedule(0.001)`: schedule used to control
strength in master. The provided function is called AFTER the incumbent
selection, but BEFORE the incumbent is replaced with candidate. The cell 
contains the candidate solution and incumbent solution from the last iteration.
See quad_scalar.jl.
- `update_incumbent_cut::Bool=true`: if true, update the incumbent cut using argmax.
"""
function sd_iteration!(cell::sdCell, scenario_list::Vector{spSmpsScenario};
     update_incumbent_cut::Bool=true, quad_scalar_schedule::Function=ConstantQuadScalarSchedule(0.1))
    # Make sure the length of the epigraph variable is the same.
    @assert(length(scenario_list) == length(cell.epivar_ref))

    # Solve the subproblem
    for i in eachindex(scenario_list)
        add_scenario!(cell.epi[i], scenario_list[i], 1.0)

        # Solve subproblem at candidate
        obj, y_opt, dual_opt = solve_problem!(cell.epi[i].prob, cell.x_candidate, scenario_list[i])
        push!(cell.dual_vertices, dual_opt)

        # Solve subproblem at incumbent
        obj, y_opt, dual_opt = solve_problem!(cell.epi[i].prob, cell.x_incumbent, scenario_list[i])
        push!(cell.dual_vertices, dual_opt)
    end

    # Remove cuts with non-zero multiplier
    if termination_status(cell.master) == OPTIMAL
        for i in eachindex(cell.epicon_ref)
            delete_index = Int[]
            for j in eachindex(cell.epicon_ref[i])
                con = cell.epicon_ref[i][j]
                if abs(dual(con)) < CUT_REMOVE_TOLERANCE
                    # Mark j for deletion
                    push!(delete_index, j)
                end
            end
            deleteat!(cell.epi[i].cuts, delete_index)
        end
    else
        # @warn("Master Problem is not solved. Skipping removing cuts.")
    end

    # Make copy of the epigraph cuts info, f_{k-1} used in incumbent selection
    # before adding any new cuts
    epi_info_last = sdEpigraphInfo[sdEpigraphInfo(epi) for epi in cell.epi]

    # Generate cuts
    for epi in cell.epi
        new_cut = build_sasa_cut(epi, cell.x_candidate, cell.dual_vertices)
        push!(epi.cuts, new_cut)
        if update_incumbent_cut
            epi.incumbent_cut = build_sasa_cut(epi, cell.x_incumbent, cell.dual_vertices)
        end
    end

    # After the cuts are generated, do we see improvement?
    # If yes, replace incumbent.
    cell.improvement_info = check_improvement(epi_info_last, cell.epi,
        cell.x_candidate, cell.x_incumbent, cell.x_ref, cell.objf_original)
    
    # Call quad_scalar_schedule to determine the next step set_optimizer
    # BEFORE the incumbent is overwritten.
    rho = quad_scalar_schedule(cell)

    if cell.improvement_info.is_improved
        cell.x_incumbent .= cell.x_candidate
    end
    
    # Solve master and store candidate solution
    add_regularization!(cell, cell.x_incumbent, rho)
    sync_cuts!(cell)
    optimize!(cell.master)
    @assert(termination_status(cell.master) == OPTIMAL)
    cell.x_candidate .= value.(cell.x_ref)

    return
end
