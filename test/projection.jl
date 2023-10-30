using JuMP, SparseArrays

struct ConstraintCoefficients{CT}
    A::Matrix{Float64}
    b::Vector{Float64}
    slacks::Vector{VariableRef}
    cons::Vector{CT}
end

struct ProjectionTemplate{CT}
    projection_model::Model
    
    # reference to the constraint d - delta == ??
    d_bind::Vector{CT}

    # constraint blocks corresponding to Ad - s = 0
    c_gt::ConstraintCoefficients
    c_eq::ConstraintCoefficients
    c_lt::ConstraintCoefficients

    # constraint blocks corresponding to l <= d <= u
    lb::Vector{Float64}
    ub::Vector{Float64}
end



# Test below

model = Model(HiGHS.Optimizer)
@variable(model, x >= 0)
@variable(model, y >= 0)
@constraint(model, mc, x + y <= 1)
@objective(model, Min, x + y)

function create_projection_template(model::Model, optimizer)
    # create a new direct model with constraints copied
    # minimize ||delta||^2
    # d - delta == given_d_bar
    # Ad - s = 0
    # s >= 0 ---- Ad >= 0
    # s == 0 ---- Ad == 0
    # s <= 0 ---- Ad <= 0
    # to disable a constraint, set the slack unsigned
    
    projection_model = direct_model(optimizer)

    N = JuMP.num_variables(model)

    # variables
    @variable(projection_model, d[1:N])
    @variable(projection_model, offset[1:N])

    # binders
    @constraint(projection_model, d_bind[i=1:N], d[i] - offset[i] == 0)

    # objective
    obj::QuadExpr = sum(offset[i]^2 for i = 1:N)
    @objective(projection_model, Min, obj)

    # create Ad -s = 0 and record coefficients objects
    c_gt = projection_constraints(model, projection_model, d, AffExpr, MOI.GreaterThan{Float64}) 
    c_eq = projection_constraints(model, projection_model, d, AffExpr, MOI.EqualTo{Float64}) 
    c_lt = projection_constraints(model, projection_model, d, AffExpr, MOI.LessThan{Float64}) 

    # record ub and lb
    lb = fill(-Inf, N)
    ub = fill(+Inf, N)
    for (i, var) in enumerate(all_variables(model))
        if has_lower_bound(var)
            lb[i] = lower_bound(var)
        end
        if has_upper_bound(var)
            ub[i] = upper_bound(var)
        end
    end

    return ProjectionTemplate(projection_model, d_bind, c_gt, c_eq, c_lt, lb, ub)
end

function projection_constraints(model::Model, projection_model::Model, d::Vector{VariableRef}, F::Type, S::Type)

    # number of constraints
    M = num_constraints(model, F, S)
    # number of variables
    N = num_variables(model)

    @assert(length(d) == N)

    # old coefficients
    A = zeros(M, N)
    b = zeros(M)

    # add slack variable
    slacks = @variable(projection_model, [i = 1:M])

    for (i, con) in enumerate(all_constraints(model, F, S))
        # record old coefficients and rhs
        for (j, var) in enumerate(all_variables(model))
            A[i, j] = normalized_coefficient(con, var)
            b[i] = normalized_rhs(con)
        end
    end

    MT = typeof(projection_model)
    # add new constraint A_i' d + s_i = 0
    new_cons = @constraint(projection_model, [i=1:M], sum(A[i,j]*d[j] for j = 1:N) + slacks[i] == 0)

    return ConstraintCoefficients(A, b, slacks, new_cons)
end

pm = create_projection_template(model, HiGHS.Optimizer())

print(pm.projection_model)
