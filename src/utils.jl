using JuMP

"""
Turn a vector into dict for reverse finding.
"""
function lookup_table(v::Vector{T})::Dict{T, Int} where T
    d = Dict{T, Int}()
    for i in eachindex(v)
        d[v[i]] = i
    end
    return d
end

"""
Create dictionary mapping from constraints to values.
"""
function value_dict(keys::Vector{VariableRef},
    vals::Vector{Float64})::Dict{VariableRef, Float64}
    @assert(length(keys) == length(vals))
    return Dict(zip(keys, vals))
end

"""
Evaluate the quadratic expression given x.
"""
function evaluate_expr(expr::Union{AffExpr, QuadExpr},
    x_ref::Vector{VariableRef}, x::Vector{Float64})::Float64
    d = value_dict(x_ref, x)
    return value(z -> d[z], expr)
end