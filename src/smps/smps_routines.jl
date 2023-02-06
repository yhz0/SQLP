"""
Change template model in place into the scenario model.
Change in rhs should be signified by col_name setting to "RHS"
Returns the changed model. If the requested row or column
does not exist, throws error.
"""
function instantiate!(sp::spStageProblem, scenario::spSmpsScenario)
    for (pos, val) in scenario
        con = constraint_by_name(sp.model, pos.row_name)
        @assert(con !== nothing, "Constraint $(pos.row_name) not in this stage problem.")
        if pos.col_name == "RHS"
            set_normalized_rhs(con, val)
        else
            var = variable_by_name(sp.model, pos.col_name)
            @assert(var !== nothing, "Variable $(pos.col_name) not in this stage problem.")
            set_normalized_coefficient(con, var, val)
        end
    end
    return
end

# """
# Instantiate with recoverable scenario.
# """
# function instantiate!(sp::spStageProblem, scenario::spSmpsScenario)
#     original_scenario::spSmpsScenario = []

#     for (pos, val) in scenario
#         con = constraint_by_name(sp.model, pos.row_name)
#         @assert(con !== nothing, "Constraint $(pos.row_name) not in this stage problem.")
#         if pos.col_name == "RHS"
#             push!(original_scenario, pos => normalized_rhs(con))
#             set_normalized_rhs(con, val)
#         else
#             var = variable_by_name(sp.model, pos.col_name)
#             push!(original_scenario, pos => normalized_coefficient(con, var))
#             @assert(var !== nothing, "Variable $(pos.col_name) not in this stage problem.")
#             set_normalized_coefficient(con, var, val)
#         end
#     end
#     return original_scenario
# end

"""
Solve the (sub)problem with given last_stage_vars. Will overwrite the coefficients in the sp,
and set fix bounds on last_stage_val.
Returns a tuple (obj, y_opt, dual_opt), which is the optimal objective, optimal (subproblem)
decision and optimal dual multipliers.
"""
function solve_problem!(sp::spStageProblem, last_stage_val::Vector{Float64}, scenario::spSmpsScenario)
    instantiate!(sp, scenario)
    fix.(sp.last_stage_vars, last_stage_val, force=true)
    optimize!(sp.model)
    if termination_status(sp.model) != OPTIMAL
        print(sp.model)
        @error "Failed to solve subproblem."
    end
    dual_opt = dual.(sp.stage_constraints)
    y_opt = value.(sp.current_stage_vars)
    obj = objective_value(sp.model)
    return obj, y_opt, dual_opt
end

"""
Evaluate a 2SP by IID sampling at a point, using N sample points.
"""
function evaluate(sp1::spStageProblem, sp2::spStageProblem, sto::spStoType, 
    x::Vector{Float64}; N::Int=10000)
    f = objective_function(sp1.model)
    s1_cost = evaluate_expr(f, sp1.current_stage_vars, x)
    s2_cost = 0.0

    # Suppress output to avoid flooding the screen.
    set_silent(sp2.model)

    for _ = 1:N
        scenario = rand(sto)
        obj, _, _ = solve_problem!(sp2, x, scenario)
        s2_cost += 1/N*obj
    end
    return s1_cost+s2_cost
end