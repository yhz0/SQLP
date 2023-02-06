const INCUMBENT_SELECTION_Q = 0.2
"""
Performs incumbent selection step.
Returns a tuple (lb_est, replace)
lb_est: lower bound estimation at current candidate
replace:
"""
function incumbent_selection(
    f_last::Vector{T1}, f_current::Vector{T2},
    x_candidate::Vector{Float64}, x_incumbent::Vector{Float64},
    x_ref::Vector{VariableRef}, common_expr::QuadExpr;
    sense::MOI.OptimizationSense = MIN_SENSE)::Tuple{Float64, Bool} where
    {T1, T2 <: Union{sdEpigraph, sdEpigraphInfo}}

    f_cand = evaluate_expr(common_expr, x_ref, x_candidate)
    f_inc = evaluate_expr(common_expr, x_ref, x_incumbent)

    lb_est = evaluate_multi_epigraph(f_current, x_candidate; sense=sense) + f_cand
    b = evaluate_multi_epigraph(f_current, x_incumbent; sense=sense) + f_inc
    cur = lb_est

    c = evaluate_multi_epigraph(f_last, x_candidate; sense=sense) + f_cand
    d = evaluate_multi_epigraph(f_last, x_incumbent; sense=sense) + f_inc
    req_improvement = INCUMBENT_SELECTION_Q * (c - d)
    req = b + req_improvement

    if sense == MIN_SENSE
        passed = cur < req
    else
        passed = cur > req
    end

    if passed
        print("+")
    else
        print("X")
    end

    # @info("Incumbent Selection: cur=$cur, req=$req, imp=$req_improvement")
    return lb_est, passed
end

const DUAL_TOLERANCE = 0.001

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
5. Generate new candidate cuts from each epigraph variable
by argmax procedure, at the candidate solution. Store it in the epigraph.
6. Do 5, update incumbent cut with incumbent solution.
7. Do incumbent selection. Replace incumbent if needed.
8. Sync the cuts to the cell's master.
Solve master to get new candidate solution.
"""
function sd_iteration!(cell::sdCell, scenario_list::Vector{spSmpsScenario}; rho::Float64 = 0.01)
    # Make sure the length of the epigraph variable is the same.
    @assert(length(scenario_list) == length(cell.epivar_ref))

    # Make copy of the epigraph cuts info, f_{k-1} used in incumbent selection
    epi_info_last = sdEpigraphInfo[sdEpigraphInfo(epi) for epi in cell.epi]

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
                if abs(dual(con)) >= DUAL_TOLERANCE
                    # Mark j for deletion
                    push!(delete_index, j)
                end
            end
            deleteat!(cell.epi[i].cuts, delete_index)
        end
    else 
        # @warn("Master Problem is not solved. Skipping removing cuts.")
    end

    # Generate cuts
    for epi in cell.epi
        new_cut = build_sasa_cut(epi, cell.x_candidate, cell.dual_vertices)
        push!(epi.cuts, new_cut)
        epi.incumbent_cut = build_sasa_cut(epi, cell.x_incumbent, cell.dual_vertices)
    end

    # After the cuts are generated, do we see improvement?
    # If yes, replace incumbent.
    lb_est, replace_incumbent = incumbent_selection(epi_info_last, cell.epi,
        cell.x_candidate, cell.x_incumbent, cell.x_ref, cell.objf_original)
    if replace_incumbent
        cell.x_incumbent .= cell.x_candidate
    end

    # Solve master
    add_regularization!(cell, cell.x_incumbent, rho)
    sync_cuts!(cell)
    optimize!(cell.master)
    
    cell.x_candidate .= value.(cell.x_ref)


    return cell.x_candidate, lb_est, replace_incumbent
end
