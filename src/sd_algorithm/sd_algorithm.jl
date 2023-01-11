using LinearAlgebra
# Implementation of Two-stage Stochastic Decomposition algorithm

"""
Internal State information during SD run.
master: master problem


"""

mutable struct sdCell
    master::Model
    
    # Reference in the order of assignment
    x_ref::Vector{VariableRef}
    root_stage_con::Vector{ConstraintRef}

    epivar_ref::Vector{VariableRef}

    sdCell() = new()
end

"""
Copy master problem and master constraints into the cell.
"""
function initialize_master!(cell::sdCell, root_prob::spStageProblem)
    cell.master, refmap = JuMP.copy_model(root_prob.model)

    cell.x_ref = [refmap[v] for v in root_prob.current_stage_vars]
    cell.root_stage_con = [refmap[c] for c in root_prob.stage_constraints]

    return
end

"""
Add epigraph variables to the master
"""
function add_epigraph_variables!(cell::sdCell, num::Int, weights::Vector{Float64})
    if sum(weights) != 1.0
        @warn("Epigraph weights does not add up to 1.0.")
    end

    @assert(length(weights) == num)

    # Register epigraph variables
    epivar = @variable(master, _eta[1:num])
    cell.epivar_ref = epivar

    # Update objective function
    f = objective_function(cell.master)
    set_objective_function(cell.master, f + dot(weights, epivar))

    return
end
