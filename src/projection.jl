using LinearAlgebra, SparseArrays, JuMP

"""
Calculate the LU factorization of [I A.T; A 0]. Returns P.
To solve the projection problem, Calculate P right div [x0; b]
"""
function factorize_projection(A::AbstractMatrix{Float64})
    m, n = size(A)
    P = factorize([UniformScaling(1.0)(n) transpose(A); A spzeros(m, m)])
    return P
end

