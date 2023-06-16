using DifferentialEquations
using StaticArrays

function kuramotoNO_static!(u, p, t)
    ω, K, N = p
    θi = reduce(hcat, repeat([u], N, 1))
    θj = transpose(θi)
    phase_diff = sin.(θj - θi)
    summed = SVector{3, Float32}(sum(phase_diff, dims = 2))
    display(typeof(ω))
    display(typeof(summed))
    SVector{3, Float32}(ω + K * summed)
end

K = [1 2 3; 4 1 6; 7 8 1]
display(K)