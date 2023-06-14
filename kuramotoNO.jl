#=
using Pkg
Pkg.add("DifferentialEquations")
Pkg.add("LinearAlgebra")
Pkg.add("Plots")
=#

using DifferentialEquations
using StaticArrays
using LinearAlgebra


function meshgrid(x, y)
    X = [x for _ in y, x in x]
    Y = [y for y in y, _ in x]
    X, Y
 end

function kuramotoNO!(du, u, p, t)
    ω, K, N = p
    θi, θj = meshgrid(u,u)
    phase_diff = sin.(θj - θi)
    summed = sum(phase_diff, dims = 2)
    for i = 1:N
        du[i] = ω[i] + K * summed[i]
    end
end

# Define initial conditions and parameters
u0 = [0.0, 1.0, 2.0]   # Initial phase values
K = 0.3e9
N = 3
p = [[6.6e9, 6.7e9, 6.2e9], K, N]    # Natural frequencies


tstart = 0.0     # Start time
tend = 1000e-9      # End time
dt = 1e-12        # Time step

tspan = (tstart, tend)  # Time span for simulation
# Define the ODE problem
prob = ODEProblem(kuramotoNO!, u0, tspan, p)

# Solve the ODE problem
sol = solve(prob, Tsit5(), abstol=1e-10,reltol=1e-10, dt=dt)

# Access the solution
t = range(tspan[1], tspan[2], length=1000)
θ1 = [sol(ti)[1] for ti in t]
θ2 = [sol(ti)[2] for ti in t]
θ3 = [sol(ti)[3] for ti in t]

f1 = diff(θ1) / (t[2] - t[1]) / 1e9
f2 = diff(θ2) / (t[2] - t[1]) / 1e9
f3 = diff(θ3) / (t[2] - t[1]) / 1e9

# Plot the results
using Plots
#=
plot(t, θ1 - θ2, label="θ")
xlabel!("Time")
ylabel!("Phase differnece")
=#
plot(f1, label="Oscillator 1")
plot!(f2, label="Oscillator 2")
plot!(f3, label="Oscillator 3")
xlabel!("Time")
ylabel!("Phase differnece")