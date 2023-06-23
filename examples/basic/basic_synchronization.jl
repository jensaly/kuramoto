using BenchmarkTools
include("../../kuramoto.jl")

N = 2
u = zeros(N)
ω = [6.6e9, 6.7e9]
K = 0.3e9
K = create_standard_K(K, N)

tstart = 0.0     # Start time
tend = 100e-9      # End time
dt = 1e-13        # Time step

tspan = (tstart, tend)  # Time span for simulation
model = Kuramoto(u, ω, K, tstart, tend, dt)

run_kuramoto(model, 1e-12, 1e-12)

plot_model_frequencies(model)