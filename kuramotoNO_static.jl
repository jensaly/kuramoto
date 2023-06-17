using DifferentialEquations
using BenchmarkTools
using StaticArrays
using Plots

 
function kuramotoNO(u, p, t)
    ω, K, N, uT, A1, A2, v1, v2 = p

    A1 = u' .- u # Phase differences
    A2 = sin.(A1)

    v1 = sum(A2, dims=2) # Summing along rows
    v2 = K * v1 # Adjancency matrix -> Interaction term

    return SVector{4, Float64}(ω .+ v2)
end

N = 4

u = @SVector zeros(N)
du = @SVector zeros(N)

ω = @SVector [6.6e9, 6.7e9, 6.2e9, 6.4e9]
K = fill(0.3e9, N, N)
K[diagind(K)] .= 0.0

K = SMatrix{N,N}(K)

u0 = @SVector zeros(N)
du = @SVector zeros(N)

# Preallocate
uT = @SVector zeros(N)
A1 = @SMatrix zeros(N, N)
A2 = @SMatrix zeros(N, N)
v1 = @SVector zeros(N)
v2 = @SVector zeros(N)

# ODE setup
tstart = 0.0     # Start time
tend = 100e-9      # End time
dt = 1e-12        # Time step

tspan = (tstart, tend)  # Time span for simulation
p = (ω, K, N, uT, A1, A2, v1, v2)

# Define the ODE problem
prob = ODEProblem(kuramotoNO, u0, tspan, p)

# Solve the ODE problem
#@btime sol = solve(prob, Tsit5(), abstol=1e-10,reltol=1e-10, dt=dt)
@btime sol = solve(prob, Tsit5(), abstol=1e-10,reltol=1e-10, dt=dt)

#=
# Access the solution
t = range(tspan[1], tspan[2], length=1000)
θ1 = [sol(ti)[1] for ti in t]
θ2 = [sol(ti)[2] for ti in t]
θ3 = [sol(ti)[3] for ti in t]
θ4 = [sol(ti)[4] for ti in t]

f1 = diff(θ1) / (t[2] - t[1]) / 1e9
f2 = diff(θ2) / (t[2] - t[1]) / 1e9
f3 = diff(θ3) / (t[2] - t[1]) / 1e9
f4 = diff(θ4) / (t[2] - t[1]) / 1e9

plot(f1, label="Oscillator 1")
plot!(f2, label="Oscillator 2")
plot!(f3, label="Oscillator 3")
plot!(f4, label="Oscillator 4")
xlabel!("Time")
ylabel!("Frequency")
=#