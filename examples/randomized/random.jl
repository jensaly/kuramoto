using BenchmarkTools
using Distributions
include("../../kuramoto.jl")

N = 100;
u = rand(N) * 2 * pi;
ω = rand(Uniform(6.2e9,7.2e9), N);
K = 0.6e9;
K = create_standard_K(K, N);

tstart = 0.0;    # Start time
tend = 1000e-9;     # End time
dt = 1e-12;        # Time step

tspan = (tstart, tend);  # Time span for simulation
model = Kuramoto(u, ω, K, tstart, tend, dt);

run_kuramoto(model, 1e-12, 1e-12)

plot_model_frequencies(model, 1, false)
