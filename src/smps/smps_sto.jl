"""
Describe a type of random variable in the SMPS sto format, INDEP section.
"""
abstract type spSmpsIndepDistribution end

"""
Scalar discrete distribution.
"""
struct spSmpsDiscreteDistribution <: spSmpsIndepDistribution
    value::Vector{Float64}
    probability::Vector{Float64}
end

"""
Scalar normal distribution with mean and variance.
"""
struct spSmpsNormalDistribution <: spSmpsIndepDistribution
    mean::Float64
    variance::Float64
end

"""
Scalar continuous uniform distribution on [a, b].
"""
struct spSmpsUniformDistribution <: spSmpsIndepDistribution
    left::Float64
    right::Float64
end

"""
Sto file representation describes distribution of coefficients.
"""
struct spStoType
    problem_name::String
    indep::Dict{spSmpsPosition, spSmpsIndepDistribution}
end

"""
Read stochastic types. There is no error checking.
"""
function read_sto(sto_path::String)::spStoType
    local lines
    open(sto_path, "r") do io
        lines = readlines(io)
    end

    local section::String
    local section_keywords::Vector{String}
    local problem_name::String
    local indep = Dict{spSmpsPosition, spSmpsIndepDistribution}()
    supported_sections::Vector{String} = ["STOCH", "INDEP", "ENDATA"]

    for line in lines
        token = split(line)

        if line[1] == ' '
            # data lines

            if section == "INDEP"
                col_name = token[1]
                row_name = token[2]
                pos = spSmpsPosition(col_name, row_name)

                # Only support univariate indep distribution and REPLACE mode for now
                if length(section_keywords) > 1
                    error("Trailing/unsupported section_keywords $section_keywords")
                end

                # use section_keywords to determine the distribution type
                if section_keywords[1] == "UNIFORM"
                    a = parse(Float64, token[3])
                    b = parse(Float64, token[4])
                    indep[pos] = spSmpsUniformDistribution(a, b)
                elseif section_keywords[1] == "NORMAL"
                    m = parse(Float64, token[3])
                    v = parse(Float64, token[4])
                    indep[pos] = spSmpsNormalDistribution(m, v)
                elseif section_keywords[1] == "DISCRETE"
                    # For discrete distribution, we first check if we have position
                    # recorded. If not, we populate it with empty distribution first.
                    if !(pos in keys(indep))
                        indep[pos] = spSmpsDiscreteDistribution(Float64[], Float64[])
                    end
                    r = Ref{spSmpsDiscreteDistribution}(indep[pos])
                    
                    v = parse(Float64, token[3])
                    p = parse(Float64, token[4])
                    push!(r[].value, v)
                    push!(r[].probability, p)

                else
                    error("Unknown or unsupported section_keywords $section_keywords")
                end
            end

        else
            #section header lines
            section = token[1]
            @assert(section in supported_sections)
            section_keywords = token[2:end]

            if section == "STOCH"
                problem_name = section_keywords[1]
            end
        end
    end

    return spStoType(problem_name, indep)
end

# Generate random variables
using Random, Distributions
function Base.rand(rng::Random.AbstractRNG, p::spSmpsDiscreteDistribution)::Float64
    dist = DiscreteNonParametric(p.value, p.probability)
    return rand(rng, dist)
end

function Base.rand(rng::Random.AbstractRNG, p::spSmpsNormalDistribution)::Float64
    dist = Normal(p.mean, sqrt(p.variance))
    return rand(rng, dist)
end

function Base.rand(rng::Random.AbstractRNG, p::spSmpsUniformDistribution)::Float64
    dist = Uniform(p.left, p.right)
    return rand(rng, dist)
end

