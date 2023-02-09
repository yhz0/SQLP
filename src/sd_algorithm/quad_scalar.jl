"""
Construct a constant learning rate schedule.
"""
function ConstantQuadScalarSchedule(reg::Float64)
    g(::sdCell)::Float64 = reg
    return g
end

"""
Construct a trust-region-like learning rate schedule.
In this schedule, we maintain :normDk_1 in extension
registry cell.ext. normDk_1 stores the norm of x_incumbent - x_candidate
from the last iteration.
"""
function AdaptiveQuadScalarSchedule(;
    min_quad_scalar::Float64=1e-3, max_quad_scalar::Float64=1e4,
    R2::Float64=0.95, R3::Float64=2.0, tolerance=1e-3
    )

    function g(cell::sdCell)::Float64
        @assert(:quad_scalar in keys(cell.ext),
        "Quad_scalar not initialized. "
        * "To use AdaptiveQuadScalarSchedule, set up cell.ext[:quad_scalar] first!")

        # Calculate ||x_inc - x_cand||^2
        normDk = 0.0
        for i in eachindex(cell.x_candidate)
            d = cell.x_incumbent[i] - cell.x_candidate[i]
            normDk += d * d
        end
        
        # Initialze normDk_1 if never done before
        if !(:normDk_1 in keys(cell.ext))
            # If there is movement, then we initialize normDk_1 to avoid
            # multiplication by zero.
            if normDk > tolerance
                cell.ext[:normDk_1] = normDk
            else
                # No movement, or uninitialized. Just use the old
                # inital quad scalar
                return cell.ext[:quad_scalar]
            end
        end
        normDk_1 = cell.ext[:normDk_1]

        if cell.improvement_info.is_improved
            # Incumbent was replaced and the step is significantly larger
            if normDk > tolerance && normDk >= R3 * normDk_1
                cell.ext[:quad_scalar] *= (R2 * R3 * normDk_1 / normDk)
            end
        else
            # Incumbent was not replaced
            cell.ext[:quad_scalar] /= R2
        end

        # In any case, clamp quad scalar to interval
        cell.ext[:quad_scalar] = min(cell.ext[:quad_scalar], max_quad_scalar)
        cell.ext[:quad_scalar] = max(cell.ext[:quad_scalar], min_quad_scalar)

        # Record normDk_1
        cell.ext[:normDk_1] = normDk

        return cell.ext[:quad_scalar]
    end

end