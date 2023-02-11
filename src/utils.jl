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

"""
Copy constraint con_ref to dest_model, with variables mapping var_map.
"""
function copy_constraint!(dest_model, con_ref, var_map)
    co = constraint_object(con_ref)
    new_con = @constraint(dest_model, map_function_vars(co.func, var_map) in co.set)
    return new_con
end

"""
Copy variables into specified model, with optionally bounds copied.
"""
function copy_variable!(dest_model, var_ref; copy_bounds::Bool=true)::VariableRef
    new_var = @variable(dest_model)
    if copy_bounds
        if has_lower_bound(var_ref)
            set_lower_bound(new_var, lower_bound(var_ref))
        end
        if has_upper_bound(var_ref)
            set_upper_bound(new_var, upper_bound(var_ref))
        end
    end
    return new_var
end

"""
Copy a function with a substitution of variables specified in
dictionary mapping var_map.
"""
function map_function_vars(expr::AffExpr, var_map)::AffExpr
    new_expr = AffExpr(expr.constant)
    for (term, val) in expr.terms
        @assert(term in keys(var_map))
        new_term = var_map[term]
        add_to_expression!(new_expr, new_term, val)
    end
    return new_expr
end
