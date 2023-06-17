using DifferentialEquations
using StaticArrays
using BenchmarkTools
using LinearAlgebra


function kuramotoNO!(du, u, p, t)
    ω, K, N, uT, A1, A2, v1, v2 = p

    uT .= u'

    A1 .= uT .- u # Finding phase differences
    A2 .= sin.(A1)

    sum!(v1, A2) # Summing along sin(phase)
    mul!(v2, K, v1) # Adjancency matrix -> Interaction term

    du .= ω .+ v2 # Step
end

N = 4

ω = rand(5.5e9:7.2e9, N)
K = fill(0.3e9, N, N)
K[diagind(K)] .= 0.0

u = zeros(N)
du = similar(u)

# Preallocate
uT = similar(u')
A1 = similar(u, N, N)
A2 = similar(u, N, N)
v1 = similar(u)
v2 = similar(u)

p = (ω, K, N, uT, A1, A2, v1, v2)

@btime kuramotoNO!($du, $u, $p, 0.0)

print("Exit code 0")