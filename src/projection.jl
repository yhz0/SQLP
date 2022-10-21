using LinearAlgebra, SparseArrays, JuMP

# """
# Given JuMP linear expression ex, return the coefficient in the order of sp.
# """
# function coefficient(ex::AffExpr, sp::spStageProblem)
#     x_coef = spzeros(length(sp.last_stage_vars))
#     y_coef = spzeros(length(sp.current_stage_vars))
    
# end

# """
# TODO: Projection onto convex set algorithm, with previous stage variable fixed.
# """
# function projection(sp::spStageProblem, prev::Vector{Float64}, x::Vector{Float64})

# end

"""
Calculate the LU factorization of [I A.T; A 0]. Returns P.
To solve the projection problem, Calculate P right div [x0; b]
"""
function factorize_projection(A::AbstractMatrix{Float64})
    m, n = size(A)
    P = factorize([UniformScaling(1.0)(n) transpose(A); A spzeros(m, m)])
    return P
end

