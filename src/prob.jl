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
