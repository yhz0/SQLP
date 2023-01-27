# Unit tests for SQLP
using Test, SQLP
# Token tests
token = SQLP._tokenize_cor(open("spInput/lands/lands.cor"))

# _parse_row_tokenss
direction, row_names = SQLP._parse_row_tokens(token["ROWS"])
@test direction == collect("NGLLLLLGGG")
@test row_names == ["OBJ", "S1C1", "S1C2", "S2C1",
"S2C2", "S2C3", "S2C4", "S2C5", "S2C6", "S2C7"]

# _parse_unique_columns
col_names = SQLP._parse_unique_columns(token["COLUMNS"])
@test col_names == [
    "X1", "X2", "X3", "X4",
    "Y11", "Y21", "Y31", "Y41",
    "Y12", "Y22", "Y32", "Y42",
    "Y13", "Y23", "Y33", "Y43"
]

# _parse_column_to_matrix
M = SQLP._parse_column_to_matrix(token["COLUMNS"], row_names, col_names)
@test count(!iszero, M) == 52

# _parse_rhs
rhs = SQLP._parse_rhs(token["RHS"], row_names)
@test rhs == Float64[0., 12, 120, 0, 0, 0, 0, 0, 3, 2]

# _parse_bounds
lb, ub = SQLP._parse_bounds(token["BOUNDS"], col_names)
@test all(lb .== 0.0)
@test all(ub .== +Inf)

## TIM section
tim = SQLP.read_tim("spInput/lands/lands.tim")
@test tim.problem_name == "LandS"
@test tim.periods[1] == SQLP.spSmpsImplicitPeriod("TIME1", SQLP.spSmpsPosition("X1", "OBJ"))
@test tim.periods[2] == SQLP.spSmpsImplicitPeriod("TIME2", SQLP.spSmpsPosition("Y11", "S2C1"))

# COR integrated tests
cor = SQLP.read_cor("spInput/lands/lands.cor")

using JuMP
# Template reading
sp1 = SQLP.get_smps_stage_template(cor, tim, 1)
@test length(all_variables(sp1.model)) == 4
@test length(all_constraints(sp1.model; include_variable_in_set_constraints=false)) == 2
sp2 = SQLP.get_smps_stage_template(cor, tim, 2)
@test length(all_variables(sp2.model)) ==  16
@test length(all_constraints(sp2.model; include_variable_in_set_constraints=false)) == 7
@test objective_function(sp2.model) != 0

# STO section
sto = SQLP.read_sto("spInput/lands/lands.sto")
@test sto.problem_name == "LandS"
pos = SQLP.spSmpsPosition("RHS", "S2C5")
@test sto.indep[pos].value == [3.0, 5.0, 7.0]
@test sto.indep[pos].probability == [0.3, 0.4, 0.3]

# Scenario generation
using Random
rng = MersenneTwister(1234)
@test rand(rng, sto.indep[pos]) in [3.0, 5.0, 7.0]

scenario = rand(sto)
@test scenario[1].second in [3.0, 5.0, 7.0]

# Test change subproblem
scenario = [SQLP.spSmpsPosition("RHS", "S2C5") => 4.0]
SQLP.instantiate!(sp2, scenario)
@test normalized_rhs(constraint_by_name(sp2.model, "S2C5")) == 4.0

# If position specified by scenario is invalid then
# it should throw an error
@test_throws AssertionError SQLP.instantiate!(sp1, scenario)
