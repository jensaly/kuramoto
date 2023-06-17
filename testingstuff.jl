using DifferentialEquations
using StaticArrays
using BenchmarkTools


u0 = SA[0.0, 1.0, 2.0, 3.0]   # Initial phase values
ω = SMatrix{4,1,Float64}([6.6e9, 6.7e9, 6.2e9, 6.4e9])
K = SA[0.0 0.3e9 0.3e9 0.3e9;
     0.3e9 0.0 0.3e9 0.3e9;
     0.3e9 0.3e9 0.0 0.3e9; 
     0.3e9 0.3e9 0.3e9 0.0]
N = length(u0)
p = SA[ω, K, N]

ω, K, N = p
display(u0')
display(u0)
interactions = K * sum(sin.(u0' .- u0))
display(typeof(ω))
display(typeof(K * summed))
display(ω + interactions)
