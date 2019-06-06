include("./setup.jl")
include("./initialize.jl")
include("./function.jl")
using .Const, .Init, .Func, LinearAlgebra

A = Array(Diagonal(fill(1.0, (5, 5))))
display(A)
println()
