module SQLP
    include("utils.jl")
    include("prob.jl")
    include("smps/smps.jl")
    include("projection.jl")

    # Algorithms
    include("sd_algorithm/sd_algorithm.jl")
end