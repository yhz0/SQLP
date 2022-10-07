# Unit tests for SQLP
using Pkg; Pkg.activate(".")

using Test, SQLP

using JuMP
model = SQLP.read_cor("lands/lands.cor")
