using SparseArrays

# Matrix coefficients in the second stage problem
struct sdSubprobCoefficients
    rhs
    transfer
    recourse
end

# Extract the r, T, W matrices from template
function extract_coefficients(prob::spStageProblem)
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
