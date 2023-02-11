using JuMP, SQLP, GLPK, Test

cor = SQLP.read_cor(joinpath("spInput", "lands", "lands.cor"))
tim = SQLP.read_tim(joinpath("spInput", "lands", "lands.tim"))
sto = SQLP.read_sto(joinpath("spInput", "lands", "lands.sto"))

# Inputs
sp1 = SQLP.get_smps_stage_template(cor, tim, 1)
sp2 = SQLP.get_smps_stage_template(cor, tim, 2)

scenarios = SQLP.spSmpsScenario[
    [SQLP.spSmpsPosition("RHS", "S2C5") => 3.0],
    [SQLP.spSmpsPosition("RHS", "S2C5") => 5.0],
    [SQLP.spSmpsPosition("RHS", "S2C5") => 7.0]
]
probs = [0.3, 0.4, 0.3]

model, annotation = SQLP.all_in_one(sp1, sp2, scenarios, probs)

# Root stage
@test length(annotation[0]) == length(sp1.current_stage_vars)
@test length(annotation[1]) == length(sp2.current_stage_vars)
@test length(keys(annotation)) == length(scenarios) + 1

# Test whether all constraints are added
@test num_constraints(model; count_variable_in_set_constraints=false) ==
    length(sp1.stage_constraints) + length(scenarios) * length(sp2.stage_constraints)

# Test whether all objectives terms are added; just test sum; they should add up same
sp1_obj_term_sum = sum(values(objective_function(sp1.model).terms))
sp2_obj_term_sum = sum(values(objective_function(sp2.model).terms))
new_obj_term_sum = sum(values(objective_function(model).terms))
@test sp1_obj_term_sum + sp2_obj_term_sum == new_obj_term_sum

set_optimizer(model, GLPK.Optimizer)
optimize!(model)
@test isapprox(objective_value(model), 381.8533333)

# # Large scale to test for annotations for CPLEX
# using CPLEX
# function add_annotation(
#     model::JuMP.Model,
#     variable_classification::Dict;
#     all_variables::Bool = true,
# )
#     num_variables = sum(length(it) for it in values(variable_classification))
#     if all_variables
#         @assert num_variables == JuMP.num_variables(model)
#     end
#     indices, annotations = CPXINT[], CPXLONG[]
#     for (key, value) in variable_classification
#         for variable_ref in value
#             push!(indices, variable_ref.index.value - 1)
#             push!(annotations, CPX_BENDERS_MASTERVALUE + key)
#         end
#     end
#     cplex = unsafe_backend(model)
#     index_p = Ref{CPXINT}()
#     CPXnewlongannotation(
#         cplex.env,
#         cplex.lp,
#         CPX_BENDERS_ANNOTATION,
#         CPX_BENDERS_MASTERVALUE,
#     )
#     CPXgetlongannotationindex(
#         cplex.env,
#         cplex.lp,
#         CPX_BENDERS_ANNOTATION,
#         index_p,
#     )
#     CPXsetlongannotations(
#         cplex.env,
#         cplex.lp,
#         index_p[],
#         CPX_ANNOTATIONOBJ_COL,
#         length(indices),
#         indices,
#         annotations,
#     )
#     return
# end

# model, annotation = SQLP.all_in_one(sp1, sp2, [rand(sto) for i in 1:10000])
# set_optimizer(model, CPLEX.Optimizer)
# MOIU.attach_optimizer(model)
# add_annotation(model, annotation)
# optimize!(model)

