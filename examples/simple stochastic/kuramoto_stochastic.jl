include("../../kuramoto.jl")
 
function kuramoto(u, p, t)
    ω, K, N, uT, A1, A2, v1, v2, nV = p

    A1 = u' .- u # Phase differences
    A2 = sin.(-A1)

    v1 = sum(A2, dims=2) # Summing along rows
    v2 = K * v1 # Adjancency matrix -> Interaction term

    return SVector{4, Float64}(ω .+ v2)
end

function kuramoto_noise(u, p, t)
    ω, K, N, uT, A1, A2, v1, v2, nV = p
    return nV # Uncorrelated white noise of power D
end

N = 4
K = 0.97e9
D = 0.000005e9

u = zeros(4)

ω = [6.6e9, 6.7e9, 6.2e9, 6.4e9]
K = create_standard_K(K, N)
#=
K[1,2] /= 3
K[2,1] /= 3
K[3,4] /= 3
K[4,3] /= 3
=#
tstart = 0.0     # Start time
tend = 100e-9      # End time
dt = 1e-12        # Time step

tspan = (tstart, tend)  # Time span for simulation
model = Kuramoto(u, ω, K, tstart, tend, dt, D)

run_kuramoto_static_stochastic(model, 1e-12, 1e-12)