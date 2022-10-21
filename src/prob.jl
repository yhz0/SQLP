"""
An SP Stage Problem or Subproblem.
"""
struct spStageProblem
    model::Model
    last_stage_vars::Vector{VariableRef}
    current_stage_vars::Vector{VariableRef}
    stage_constraints::Vector{ConstraintRef}
end
