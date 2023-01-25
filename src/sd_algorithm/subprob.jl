using SparseArrays

# Matrix coefficients in the second stage problem
mutable struct sdSubprobCoefficients
    rhs
    transfer
    recourse

    # Lookup table from row or column names to indices
    col_lookup
    row_lookup
end

# Extract the r, T, W matrices from template
function extract_coefficients(prob::spStageProblem)::sdSubprobCoefficients

    # If the problem is non-linear or contains non-trivial bounds
    # in second stages then it should not be supported.
    for var in prob.current_stage_vars
        if has_upper_bound(var) && upper_bound(var) != 0.0
            @warn "$var has non-trivial upper bound."
        end
        if has_lower_bound(var) && lower_bound(var) != 0.0
            @warn "$var has non-trivial lower bound."
        end
    end

    for (F, S) in list_of_constraint_types(prob.model)
        if F != AffExpr && F != VariableRef
            @warn "Second stage problem contains nonlinear constraint type $F"
        end
    end

    T = spzeros(length(prob.stage_constraints), length(prob.last_stage_vars))
    for (i, con) in enumerate(prob.stage_constraints)
        for (j, var) in enumerate(prob.last_stage_vars)
            nc = normalized_coefficient(con, var)
            if !isapprox(nc, 0.0)
                T[i, j] = nc
            end
        end
    end

    r = spzeros(length(prob.stage_constraints))
    for (i, con) in enumerate(prob.stage_constraints)
        nc = normalized_rhs(con)
        if !isapprox(nc, 0.0)
            r[i] = nc
        end
    end

    W = spzeros(length(prob.stage_constraints), length(prob.current_stage_vars))
    for (i, con) in enumerate(prob.stage_constraints)
        for (j, var) in enumerate(prob.current_stage_vars)
            nc = normalized_coefficient(con, var)
            if !isapprox(nc, 0.0)
                W[i, j] = nc
            end
        end
    end

    col_names = name.(prob.last_stage_vars)
    col_lookup = lookup_table(col_names)

    row_names = name.(prob.stage_constraints)
    row_lookup = lookup_table(row_names)

    return sdSubprobCoefficients(r, T, W, col_lookup, row_lookup)
end

"""
Modify the coefficients according to the scenario.
"""
function modify_coefficients!(coef::sdSubprobCoefficients, smps_scenario)
    for i in eachindex(smps_scenario)
        pos::spSmpsPosition = smps_scenario[i].first
        val::Float64 = smps_scenario[i].second
        
        row_index = coef.row_lookup[pos.row_name]
        if pos.col_name == "RHS" || pos.col_name == "rhs"
            coef.rhs[row_index] = val
        else
            col_index = coef.col_lookup[pos.col_name]
            coef.transfer[row_index, col_index] = val
        end
    end

    return coef
end

"""
Evaluate the dual solution dual at x, with the scenario set to omega.
This assumes that all bounds-related dual variables are trivial.
"""
function eval_dual(coef::sdSubprobCoefficients, x, dual, smps_scenario)::Float64
    # set the coefficient to omega first
    modify_coefficients!(coef, smps_scenario)
    return dot(dual, coef.rhs - coef.transfer * x)
end