"""
Turn a vector into dict for reverse finding.
"""
function lookup_table(v::Vector{T})::Dict{T, Int} where T
    d = Dict{T, Int}()
    for i in eachindex(v)
        d[v[i]] = i
    end
    return d
end
