
"""
Change template model in place into the scenario model.
Change in rhs should be signified by col_name setting to "RHS"
Returns the changed model. If the requested row or column
does not exist, throws error.
Returns a fake scenario that can be used to recover the original sp.
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
    return sp
end

# """
# Solve the subproblem with given last_stage_vars. Will not overwrite the coefficients in the sp.
# """
# function solve_subproblem(sp::spStageProblem, last_stage_val::Vector{Float64}, scenario::spSmpsScenario)
    
# end