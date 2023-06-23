include("../../kuramoto.jl")

N = 2
K = 3.0e9
D = 0.00000001e9

display(D)

u = zeros(N)

ω = [6.5e9, 8.8e9]
K = create_standard_K(K, N)

tstart = 0.0     # Start time
tend = 10e-9      # End time
dt = 1e-16        # Time step

tspan = (tstart, tend)  # Time span for simulation
model = Kuramoto(u, ω, K, tstart, tend, dt, D)

run_kuramoto_static_stochastic(model, 1e-12, 1e-12)

plot_model_frequencies(model, 200)

model2 = Kuramoto(u, ω, K, tstart, tend, dt, 0)

run_kuramoto_static(model2, 1e-12, 1e-12)

plot_model_frequencies(model2, 10)