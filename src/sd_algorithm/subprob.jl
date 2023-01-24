using SparseArrays

# Matrix coefficients in the second stage problem
struct sdSubprobCoefficients
    rhs
    transfer
    recourse
end

# Extract the r, T, W matrices from template
function extract_coefficients(prob::spStageProblem)

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

    return sdSubprobCoefficients(r, T, W)
end
