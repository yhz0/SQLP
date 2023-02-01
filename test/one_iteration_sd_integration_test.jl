using SQLP, MathOptInterface, JuMP, CPLEX, Test

optimizer = CPLEX.Optimizer
# Load files
cor = SQLP.read_cor(joinpath("spInput", "lands", "lands.cor"))
tim = SQLP.read_tim(joinpath("spInput", "lands", "lands.tim"))
sto = SQLP.read_sto(joinpath("spInput", "lands", "lands.sto"))

sp1 = SQLP.get_smps_stage_template(cor, tim, 1)
sp2 = SQLP.get_smps_stage_template(cor, tim, 2)

# Set up cell
cell = SQLP.sdCell(sp1)
epi1 = SQLP.sdEpigraph(sp2, 0.5, 0.0)
epi2 = SQLP.sdEpigraph(sp2, 0.5, 0.0)
SQLP.bind_epigraph!(cell, epi1)
SQLP.bind_epigraph!(cell, epi2)

set_optimizer(cell.master, optimizer)
set_silent(cell.master)
for epi in cell.epi
    set_optimizer(epi.prob.model, optimizer)
    set_silent(epi.prob.model)
end


# A starting solution
x0 = [3., 3, 3, 3]
@test SQLP.check_first_stage_feasible(sp1, x0; optimizer)
cell.x_incumbent .= x0
cell.x_candidate .= x0

# Populate with initial samples
for i = 1:1000
    SQLP.add_scenario!(cell.epi[1], rand(sto))
    SQLP.add_scenario!(cell.epi[2], rand(sto))
end

# End of setup. Start SD procedure.

# Test manual "SD" iteration
begin 
    # Make copy of the epigraph cuts info, f_{k-1} used in incumbent selection
    epi_info_last = [SQLP.sdEpigraphInfo(epi) for epi in cell.epi]

    scenario_list = [rand(sto), rand(sto)]

    # Solve the subproblem
    for i in eachindex(scenario_list)
        SQLP.add_scenario!(cell.epi[i], scenario_list[i], 1.0)

        # Solve subproblem at candidate
        obj, y_opt, dual_opt = SQLP.solve_problem!(cell.epi[i].prob, cell.x_candidate, scenario_list[i])
        push!(cell.dual_vertices, dual_opt)

        # Test that argmax gives optimal objective
        delta = SQLP.delta_coefficients(cell.epi[i].subproblem_coef, scenario_list[i])
        val, _ = SQLP.argmax_procedure(cell.epi[i].subproblem_coef,
            [delta], cell.x_candidate, cell.dual_vertices)
        @test val[1] ≈ obj
        
        # Solve subproblem at incumbent
        obj, y_opt, dual_opt = SQLP.solve_problem!(cell.epi[i].prob, cell.x_incumbent, scenario_list[i])
        push!(cell.dual_vertices, dual_opt)

        # Test at incumbent
        val, _ = SQLP.argmax_procedure(cell.epi[1].subproblem_coef,
            [delta], cell.x_incumbent, cell.dual_vertices)
        @test val[1] ≈ obj
    end

    # Generate cuts
    for epi in cell.epi
        new_cut = SQLP.build_sasa_cut(epi, cell.x_candidate, cell.dual_vertices)
        push!(epi.cuts, new_cut)
        epi.incumbent_cut = SQLP.build_sasa_cut(epi, cell.x_incumbent, cell.dual_vertices)
    end

    # Test epigraphs evaluated correctly
    v1 = SQLP.evaluate_epigraph(epi1, cell.x_candidate)
    v2 = SQLP.evaluate_epigraph(epi2, cell.x_candidate)
    v3 = SQLP.evaluate_multi_epigraph(cell.epi, cell.x_candidate)
    @test v1 + v2 ≈ v3

    # print(cell.master)
    # Solve master
    rho = 10.0
    SQLP.add_regularization!(cell, cell.x_incumbent, rho)
    SQLP.sync_cuts!(cell)
    optimize!(cell.master)

    @show value.(cell.x_ref) - cell.x_incumbent
    cell.x_candidate .= value.(cell.x_ref)

    # # Incumbent selection
    replace_incumbent = SQLP.incumbent_selection(epi_info_last, cell.epi,
        cell.x_candidate, cell.x_incumbent)
    if replace_incumbent
        cell.x_incumbent .= cell.x_candidate
    end

    epi_info_last = [SQLP.sdEpigraphInfo(epi) for epi in cell.epi]

end

