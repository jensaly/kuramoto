using DifferentialEquations
using StaticArrays
using BenchmarkTools

u = SA[0,1,2,3]

@btime stack([u for j=1:4])

@btime reduce(hcat, repeat([u], N, 1))