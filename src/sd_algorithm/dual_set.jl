"""
Number of significant digits to compare for a dual vertex.
"""
const SIGNIFICANT_DIGITS::Int = 16

"""
Representation of a dual vertex in D dimension.
"""
struct sdDualVertex
    hash::UInt
    data::Vector{Float64}
end

# For overloading
import Base.hash, Base.isequal, Base.push!, Base.show, Base.iterate

"""
Recursive hashing for dual vertex type.
"""
function Base.hash(dv::sdDualVertex, h::UInt)::UInt
    return h + dv.hash
end

function Base.isequal(dv1::sdDualVertex, dv2::sdDualVertex)::Bool
    # Length and hash check
    if length(dv1.data) != length(dv2.data) || dv1.hash != dv2.hash
        return false
    end

    # approximate equal check
    for i in eachindex(dv1.data)
        r1 = round(dv1.data[i]; base=2, sigdigits=SIGNIFICANT_DIGITS)
        r2 = round(dv2.data[i]; base=2, sigdigits=SIGNIFICANT_DIGITS)
        if r1 != r2
            return false
        end
    end

    return true
end

"""
Calculates one norm of the vector, and then round the sum
to get the hash of the vector.
"""
function hash_dual_vector(vec::Vector{Float64})::UInt
    mysum::Float64 = 0.0
    for v in vec
        mysum += abs(v)
    end
    h = reinterpret(UInt, round(mysum; base=2, sigdigits=SIGNIFICANT_DIGITS))
    return h
end

"""
Convenient method for constructing sdDualVertex.
This function copies reference to vec as is. In-place
modification to vec will cause undefined behavior.
"""
function sdDualVertex(vec::Vector{Float64})
    h = hash_dual_vector(vec)
    return sdDualVertex(h, vec)
end

"""
Container used to store dual vertices.
Vectors are hashed via truncated representation in their one norm.
"""
struct sdDualVertexSet
    data::Vector{sdDualVertex}
end

"""
Construct an empty dual vertex set.
"""
function sdDualVertexSet() 
    return sdDualVertexSet([])
end

"""
Add a copy of dual vertex to the set. The vector is hashed with one norm.
Return true if successful, or false if the vec already exists.
"""
function Base.push!(dvs::sdDualVertexSet, new_vec::Vector{Float64})
    temp = sdDualVertex(new_vec)
    for v in dvs.data
        if isequal(temp, v)
            return dvs
        end
    end
    push!(dvs.data, temp)
    return dvs
end

"""
Construct dual vertex set with given existing dual extreme points.
"""
function sdDualVertexSet(data::Union{Vector{Vector{Float64}}, Set{Vector{Float64}}})
    dvs = sdDualVertexSet()
    for d::Vector{Float64} in data
        push!(dvs, d)
    end
    return dvs
end

"""
Size of the container.
"""
function Base.length(dvs::sdDualVertexSet)
    return length(dvs.data)
end

"""
Iterate utilities to support for loops.
"""
function Base.iterate(dvs::sdDualVertexSet, state=1)
    if state > length(dvs)
        return nothing
    else
        return (dvs.data[state].data, state+1)
    end
end

"""
Type hint for sdDualVertexSet
"""
Base.eltype(::Type{sdDualVertexSet}) = Vector{Float64}
