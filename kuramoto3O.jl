using DifferentialEquations

function kuramoto3O!(du, u, p, t)
    θ1, θ2, θ3 = u
    ω1, ω2, ω3, K, N = p
    
    du[1] = ω1 + K * (sin(θ2 - θ1) + sin(θ3 - θ1))
    du[2] = ω2 + K * (sin(θ1 - θ2) + sin(θ3 - θ2))
    du[3] = ω3 + K * (sin(θ1 - θ3) + sin(θ2 - θ3))
end

# Define initial conditions and parameters
u0 = [0.0, 0.0, 0.0]   # Initial phase values
K = 0.1e9
N = 3
p = [6.6e9, 6.7e9, 0e9, K, N]    # Natural frequencies

tstart = 0.0     # Start time
tend = 1000e-9      # End time
dt = 1e-12        # Time step

tspan = (tstart, tend)  # Time span for simulation
# Define the ODE problem
prob = ODEProblem(kuramoto3O!, u0, tspan, p)

# Solve the ODE problem
sol = solve(prob, Tsit5(), abstol=1e-10,reltol=1e-10, dt=dt)

# Access the solution
t = range(tspan[1], tspan[2], length=10000)
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
#plot!(f3, label="Oscillator 3")
xlabel!("Time")
ylabel!("Phase differnece")