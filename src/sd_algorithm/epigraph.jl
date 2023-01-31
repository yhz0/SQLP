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
    prob::spStageProblem

    # Subproblem template coefficients shared by all scenarios.
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
    # TODO: change order of data structure member
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
"""
function build_sasa_cut(epi::sdEpigraph, x::Vector{Float64},
    dual_vertices::Union{Set{Vector{Float64}}, Vector{Vector{Float64}}})::sdCut
    # TODO
    error("Unimplemented")
end

"""
Evaluate the pointwise min(or max) at a certain x, given discount, and lower bound.
"""
function evaluate_epigraph(epi::sdEpigraph, x::Vector{Float64})
    sense = objective_sense(epi.prob.model)
    @assert(sense == MIN_SENSE || sense == MAX_SENSE)
    
    best_val = epi.lower_bound

    for cut::sdCut in epi.cuts
        discount = cut.weight_mark / epi.total_scenario_weight
        val = discount * (cut.alpha + dot(cut.beta, x)) + (1-discount) * epi.lower_bound
        if sense == MIN_SENSE && val > best_val
            best_val = val
        elseif sense == MAX_SENSE && val < best_val
            best_val = val
        end
    end

    # Eval at Incumbent cut if needed
    if epi.incumbent_cut !== nothing
        val = epi.incumbent_cut.alpha + dot(epi.incumbent_cut.beta, x)
        if sense == MIN_SENSE && val > best_val
            best_val = val
        elseif sense == MAX_SENSE && val < best_val
            best_val = val
        end
    end

    return best_val
end
