"""
A representation of a cut under an epigraph in SD.
Stands for η≥α+βx. Note that cuts stored here are never scaled.
"""
struct sdCut
    # Coefficients stands for alpha + beta' x
    alpha::Float64
    beta::Vector{Float64}
    
    # Records when this cut is generated to scale the cuts.
    weight_mark::Float64
end

"""
A convex piecewise function approximation.
"""
mutable struct sdEpigraph
    # Subproblem corresponding to this epigraph variable.
    # The model under this stage problem should be unique to this instance.
    # and may be changed externally (e.g. solve_subproblem)
    prob::spStageProblem

    # Subproblem template coefficients shared by all scenarios.
    # Once Initialized this coefficient is guaranteed not to change.
    subproblem_coef::sdSubprobCoefficients

    # Epigraph weights relative to the master problem.
    objective_weight::Float64

    # Absolute lower bound of the second stage
    lower_bound::Float64

    # Scenario weights
    scenario_weight::Vector{Float64}
    # Total scenario weights
    total_scenario_weight::Float64

    # Scenarios
    scenario_list::Vector{spSmpsScenario}
    scenario_delta::Vector{sdDeltaCoefficients}

    # Store cuts
    cuts::Vector{sdCut}
    incumbent_cut::Union{Nothing, sdCut}
end

"""
Create a new epigraph with the specified weight and subproblem template.
Make sure the subprob is not reused in other epigraph variables
so coefficients change does not overwrite each other in parallel version.
"""
function sdEpigraph(
    prob::spStageProblem,
    objective_weight::Float64,
    lower_bound::Float64)
    
    new_prob = copy(prob)
    coef = extract_coefficients(new_prob)
    return sdEpigraph(new_prob, coef, objective_weight, lower_bound,
        [], 0.0, [], [], [], nothing)
end

"""
Print an epigraph variable.
"""
function Base.show(io::IO, epi::sdEpigraph)
    objw = epi.objective_weight
    lb = epi.lower_bound
    sc_cnt = length(epi.scenario_list)
    sc_w = epi.total_scenario_weight
    cut_cnt = length(epi.cuts)
    inc_cut = (epi.incumbent_cut !== nothing)
    print(io, "sdEpigraph objw=$objw lb=$lb sc_cnt=$sc_cnt sc_w=$sc_w cut_cnt=$cut_cnt inc_cut=$inc_cut")
    return
end


"""
Add a scenario to epigraph variable epi.
"""
function add_scenario!(
    epi::sdEpigraph,
    scenario::spSmpsScenario,
    weight::Float64 = 1.0)

    # Record raw scenario
    push!(epi.scenario_list, scenario)
    push!(epi.scenario_weight, weight)
    epi.total_scenario_weight += weight

    # Calculate delta coefficients, store.
    delta = delta_coefficients(epi.subproblem_coef, scenario)
    push!(epi.scenario_delta, delta)
    
    return
end

"""
Add one cut to given master problem. Returns reference to the constraint added.
"""
function add_cut_to_master!(master::Model, cut::sdCut, 
    epi_ref::VariableRef, x_ref::Vector{VariableRef},
    discount::Float64, lower_bound::Float64)::ConstraintRef
    
    new_alpha = discount * cut.alpha + (1-discount) * lower_bound
    new_beta = discount * cut.beta

    if objective_sense(master) == MOI.MIN_SENSE
        con = @constraint(master, epi_ref >= new_alpha + dot(new_beta, x_ref))
    elseif objective_sense(master) == MOI.MAX_SENSE
        con = @constraint(master, epi_ref <= new_alpha + dot(new_beta, x_ref))
    else
        error("Unknown master sense. Master should either be MIN or MAX problem.")
    end
        
    return con
end

"""
Build a cut at given last_stage_var x and partial dual_vertices by argmax procedure.
Calculates α+βx=sum(pj π(x, ωj)(rj - Tj x)).
pi is the relatative weight of a sample in this epigraph, proportional
to weight_i and sums to 1.
"""
function build_sasa_cut(epi::sdEpigraph, x::Vector{Float64},
    dual_vertices::sdDualSet)::sdCut
    max_val, max_arg = argmax_procedure(epi.subproblem_coef, epi.scenario_delta,
        x, dual_vertices; sense=objective_sense(epi.prob.model))

    alpha::Float64 = 0.0
    beta::Vector{Float64} = zeros(length(x))
    val::Float64 = 0.0

    for i in eachindex(dual_vertices)
        # dual vertex chosen
        dual = max_arg[i]
        # weight
        p = epi.scenario_weight[i] / epi.total_scenario_weight

        alpha += p * dot(dual[], epi.subproblem_coef.rhs + epi.scenario_delta[i].delta_rhs)
        beta += p * (epi.subproblem_coef.transfer + epi.scenario_delta[i].delta_transfer)' * dual[]
        val += p * max_val[i]
    end

    return sdCut(alpha, beta, epi.total_scenario_weight)
end

"""
Structure that stores all the information to recover the
convex piecewise function approximation.
"""
struct sdEpigraphInfo
    objective_weight::Float64
    cuts::Vector{sdCut}
    incumbent_cut::Union{sdCut, Nothing}
    total_scenario_weight::Float64
    lower_bound::Float64
end

"""
Extract only the necessary info from the epigraph variables.
"""
function sdEpigraphInfo(epi::sdEpigraph)
    return sdEpigraphInfo(
        epi.objective_weight,
        copy(epi.cuts),
        epi.incumbent_cut,
        epi.total_scenario_weight,
        epi.lower_bound
    )
end

"""
Evaluate the pointwise max(or min) at a certain x, and lower bound.
This does not apply the weight of the epigraph.
"""
function evaluate_epigraph(cuts::Vector{sdCut}, incumbent_cut::Union{sdCut, Nothing},
    x::Vector{Float64}, total_scenario_weight::Float64, lower_bound::Float64;
    sense::MOI.OptimizationSense=MIN_SENSE)
    best_val = lower_bound

    for cut::sdCut in cuts
        discount = cut.weight_mark / total_scenario_weight
        val = discount * (cut.alpha + dot(cut.beta, x)) + (1-discount) * lower_bound
        if sense == MIN_SENSE && val > best_val
            best_val = val
        elseif sense == MAX_SENSE && val < best_val
            best_val = val
        end
    end

    # Eval at Incumbent cut if needed
    if incumbent_cut !== nothing
        val = incumbent_cut.alpha + dot(incumbent_cut.beta, x)
        if sense == MIN_SENSE && val > best_val
            best_val = val
        elseif sense == MAX_SENSE && val < best_val
            best_val = val
        end
    end

    return best_val
end

"""
Evaluate the epigraph given the info. This includes the weight of the epigraph.
"""
function evaluate_epigraph(epi_info::sdEpigraphInfo, x::Vector{Float64};
    sense::MOI.OptimizationSense)
    return epi_info.objective_weight * evaluate_epigraph(epi_info.cuts, epi_info.incumbent_cut,
    x, epi_info.total_scenario_weight, epi_info.lower_bound; sense=sense)
end

"""
Shortcut for evaluating the epigraph at the info level.
"""
function evaluate_epigraph(epi::sdEpigraph, x::Vector{Float64};
    sense::MOI.OptimizationSense=MIN_SENSE)
    return evaluate_epigraph(sdEpigraphInfo(epi), x; sense=sense)
end

"""
Evaluate weighted sum of epigraph at given x.
"""
function evaluate_multi_epigraph(v_epi::Vector{T}, x::Vector{Float64};
    sense::MOI.OptimizationSense=MIN_SENSE) where T <: Union{sdEpigraph, sdEpigraphInfo}
    return sum(evaluate_epigraph(epi, x; sense=sense) for epi in v_epi)
end