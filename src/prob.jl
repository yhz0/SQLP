using JuMP

"""
An SP Stage Problem or Subproblem.
model: a JuMP model that refers to the current stage.
last_stage_vars: reference to variables in the last stage. Empty for root stage.
current_stage_vars: reference to variables in the current stage.
stage_constraints: reference to constraints in the current stage.
"""
struct spStageProblem
    model::Model
    last_stage_vars::Vector{VariableRef}
    current_stage_vars::Vector{VariableRef}
    stage_constraints::Vector{ConstraintRef}
end

"""
Check feasibility of a first stage decision.
"""
function check_first_stage_feasible(sp1::spStageProblem, x::Vector{Float64}; optimizer)
    model, refmap = copy_model(sp1.model)
    fix.(refmap[sp1.current_stage_vars], x, force=true)
    set_objective_sense(model, FEASIBILITY_SENSE)
    set_optimizer(model, optimizer)
    set_silent(model)
    optimize!(model)
    if termination_status(model) == OPTIMAL
        return true
    else
        return false
    end
end

"""
Copy the problem and its underlying model.
"""
function Base.copy(prob::spStageProblem)::spStageProblem
    new_model, refmap = copy_model(prob.model)
    last_stage_vars = [refmap[v] for v in prob.last_stage_vars]
    current_stage_vars = [refmap[v] for v in prob.current_stage_vars]
    stage_constraints = [refmap[v] for v in prob.stage_constraints]
    return spStageProblem(
        new_model, last_stage_vars, current_stage_vars, stage_constraints)
end