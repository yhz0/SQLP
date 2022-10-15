using JuMP

"""
An SP problem defined by smps format.
"""
struct spSmpsProblem <: spProblem
    cor::spCorType
    tim::spTimType
end

"""
A realized scenario
"""
abstract type spNoise end

"""
Get a template of the matrix.
"""
function get_stage_template(prob::spSmpsProblem, stage::Int, in_state)::Model
    cor = prob.cor
    model = Model()
    vars = VariableRef[]
    cons = ConstraintRef[]
    
    # Add variables
    for i in eachindex(cor.col_names)
        var::VariableRef = @variable(model)
        set_name(var, cor.col_names[i])
        set_lower_bound(var, cor.lower_bound[i])
        set_upper_bound(var, cor.upper_bound[i])
        push!(vars, var)
    end

    col_name_mapping = get_name_mapping(cor.col_names)
    # Add constraints
    for i in eachindex(cor.row_names)
        
    end
end
