const INCUMBENT_SELECTION_Q = 0.2

"""
Return result for incumbent selection information.
"""
struct sdImprovementInfo
    candidate_estimation::Float64
    incumbent_estimation::Float64
    required_improvement::Float64
    is_improved::Bool
end

"""
Performs incumbent selection step.
Returns a tuple (lb_est, replace)
lb_est: lower bound estimation at current candidate
replace:
"""
function check_improvement(
    f_last::Vector{T1}, f_current::Vector{T2},
    x_candidate::Vector{Float64}, x_incumbent::Vector{Float64},
    x_ref::Vector{VariableRef}, common_expr::QuadExpr;
    sense::MOI.OptimizationSense = MIN_SENSE)::sdImprovementInfo where
    {T1, T2 <: Union{sdEpigraph, sdEpigraphInfo}}

    f_cand = evaluate_expr(common_expr, x_ref, x_candidate)
    f_inc = evaluate_expr(common_expr, x_ref, x_incumbent)

    candidate_estimation = evaluate_multi_epigraph(f_current, x_candidate; sense=sense) + f_cand
    incumbent_estimation = evaluate_multi_epigraph(f_current, x_incumbent; sense=sense) + f_inc

    last_candidate_estimation = evaluate_multi_epigraph(f_last, x_candidate; sense=sense) + f_cand
    last_incumbent_estimation = evaluate_multi_epigraph(f_last, x_incumbent; sense=sense) + f_inc

    required_improvement = INCUMBENT_SELECTION_Q * (last_candidate_estimation - last_incumbent_estimation)
    req = incumbent_estimation + required_improvement

    if sense == MIN_SENSE
        is_improved = candidate_estimation < req
    else
        is_improved = candidate_estimation > req
    end

    info = sdIncumbentSelectionInfo(
        candidate_estimation, incumbent_estimation,
        required_improvement, is_improved
    )
    return info
end
