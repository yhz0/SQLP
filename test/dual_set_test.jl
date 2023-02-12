using Test, SQLP
v1 = [1., 2, 3]
v2 = [1.0000000001, 2, 3]
v3 = [4., 5, 6]
v4 = [4., 5, 6, 7]
v5 = [3., 2, 1]

# Equality test
@test isequal(SQLP.sdDualVertex(v1), SQLP.sdDualVertex(v2))
@test isequal(SQLP.sdDualVertex(v3), SQLP.sdDualVertex(v3))
@test !isequal(SQLP.sdDualVertex(v1), SQLP.sdDualVertex(v3))
@test !isequal(SQLP.sdDualVertex(v3), SQLP.sdDualVertex(v4))
@test !isequal(SQLP.sdDualVertex(v5), SQLP.sdDualVertex(v1))

# Remove duplication
dvs = SQLP.sdDualVertexSet()
push!(dvs, v1)
@test length(dvs) == 1
push!(dvs, v2)
@test length(dvs) == 1
push!(dvs, v3)
@test length(dvs) == 2
push!(dvs, v4)
@test length(dvs) == 3
push!(dvs, v5)
@test length(dvs) == 4

# Initialize with De-duplication
dvs = SQLP.sdDualVertexSet([v1,v2,v3,v4,v5])
@test length(dvs) == 4

# Iterate
@test length(collect(dvs)) == 4

# for dv in dvs
#     @show dv
# end
