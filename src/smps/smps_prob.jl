using JuMP

"""
An SP problem defined by smps format.
"""
struct spSmpsProblem <: spProblem
    cor::spCorType
    tim::spTimType
end

"""
An SP Stage Problem or Subproblem.
"""
struct spStageProblem
    model::Model
    last_stage_vars::Union{Vector{VariableRef}, Nothing}
    current_stage_vars::Vector{VariableRef}
    stage_constraints::Vector{ConstraintRef}
end

"""
Get a template model from SMPS cor and tim partitions, returns spStageProblem template.
"""
function get_smps_stage_template(cor::spCorType, tim::spTimType, stage::Int)::spStageProblem
    model = Model()

    @assert(1 <= stage <= length(tim.periods))

    # Determine starting column, to point to stage-1 problem
    start_col_index = if stage == 1
        1
    else
        cor.col_mapping[tim.periods[stage-1].position.col_name]
    end

    # Determine ending column.
    # If terminal stage then add until end
    # otherwise until one less than next stage
    end_col_index = if stage < length(tim.periods)
        cor.col_mapping[tim.periods[stage+1].position.col_name] - 1
    else 
        length(cor.col_names)
    end

    current_stage_start_col_index = cor.col_mapping[tim.periods[stage].position.col_name]

    last_stage_vars = VariableRef[]
    current_stage_vars = VariableRef[]

    # Add variables from stage-1 and this stage
    for i = start_col_index:end_col_index
        var::VariableRef = @variable(model)
        set_name(var, cor.col_names[i])
        if cor.upper_bound[i] != Inf
            set_upper_bound(var, cor.upper_bound[i])
        end
        if cor.lower_bound[i] != -Inf
            set_lower_bound(var, cor.lower_bound[i])
        end

        if i >= current_stage_start_col_index
            push!(current_stage_vars, var)
        else
            push!(last_stage_vars, var)
        end
    end

    # We do things on the first objective row, and only add the current stage
    obj = cor.template_matrix[1, current_stage_start_col_index:end_col_index]' * current_stage_vars
    set_objective(model, MOI.MIN_SENSE, obj)

    # We checked that the first row is the objective row, so we will disregard that now.
    # Now copy all other rows
    start_row_index = if stage == 1
        2
    else
        cor.row_mapping[tim.periods[stage].position.row_name]
    end

    # Determine ending column.
    # If terminal stage then add until end
    # otherwise until one less than next stage
    end_row_index = if stage < length(tim.periods)
        cor.row_mapping[tim.periods[stage+1].position.row_name] - 1
    else 
        length(cor.row_names)
    end

    cons = ConstraintRef[]
    for i in start_row_index:end_row_index
        func = (
            cor.template_matrix[i, start_col_index:end_col_index]'
            * vcat(last_stage_vars, current_stage_vars)
        )
        rhs = cor.rhs[i]
        
        dir = cor.directions[i]
        if dir == 'G'
            con = @constraint(model, func >= rhs)
        elseif dir == 'L'
            con = @constraint(model, func <= rhs)
        elseif dir == 'E'
            con = @constraint(model, func == rhs)
        else
            error("Unknown direction $dir")
        end
        set_name(con, cor.row_names[i])
        push!(cons, con)
    end

    return spStageProblem(model, last_stage_vars, current_stage_vars, cons)
end

"""
Change template model in place into the scenario model.
Change in rhs should be signified by col_name setting to "RHS"
Returns the changed model. If the requested row or column
does not exist, throws error.
"""
function instantiate!(sp::spStageProblem, scenario::Vector{Pair{spSmpsPosition, Float64}})::spStageProblem
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
