using JuMP

"""
An SP Stage Problem or Subproblem.
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
    optimize!(model)
    if termination_status(model) == OPTIMAL
        return true
    else
        return false
    end
end