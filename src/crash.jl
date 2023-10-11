"""
These are heuristic methods that are used to produce a solution that is "acceptable"
for a starting point.
"""

using JuMP


"""
Formulate all-in-one problem with root stage sp1, second stage sp2,
scenarios and specified probabilities.
Returns the generated model and variable annotations,
possibly used for Benders. The variable annotation is a dictionary of form:
    0 => [x1, x2, ...] for root stage
    1 => [y11, y12, y13, ...] for second stage, scenario 1
    2 => [y21, y22, y23, ...] for second stage, scenario 2 and so on. 
"""
function all_in_one(sp1::spStageProblem, sp2::spStageProblem,
    scenarios::Vector{spSmpsScenario}, probs::Union{Nothing, Vector{Float64}} = nothing)

    # Some sanity checks
    @assert(length(sp2.last_stage_vars) == length(sp1.current_stage_vars))

    # Copy the root stage first
    model, main_map = copy_model(sp1.model)
    root_stage_vars = getindex.(main_map, sp1.current_stage_vars)

    # Add root stage annotations
    annotations = Dict{Int, Vector{VariableRef}}(0 => root_stage_vars)

    # Record root stage objective
    objf = objective_function(model)

    for s in eachindex(scenarios)
        if probs === nothing
            prob = 1.0/length(scenarios)
        else
            prob = probs[s]
        end
        scenario = scenarios[s]

        TwoSD.instantiate!(sp2, scenario)

        # Create variables for second stage; copy bounds
        second_stage_vars = copy_variable!.(model, sp2.current_stage_vars; copy_bounds=true)
        # Add to annotations
        push!(annotations, s => second_stage_vars)

        # Create variable
        secondary_map = Dict{VariableRef, VariableRef}()
        for i in eachindex(sp2.last_stage_vars)
            secondary_map[sp2.last_stage_vars[i]] = root_stage_vars[i]
        end
        for i in eachindex(second_stage_vars)
            secondary_map[sp2.current_stage_vars[i]] = second_stage_vars[i]
        end

        # Accumulate Objective with Weights
        oterms = map_function_vars(objective_function(sp2.model), secondary_map)
        add_to_expression!(objf, prob, oterms)
        
        # Copy Constraints from this scenario
        for con in all_constraints(sp2.model; include_variable_in_set_constraints=false)
            copy_constraint!(model, con, secondary_map)
        end

    end

    # Set objective with accumulated function
    set_objective_function(model, objf)

    return model, annotations
end
