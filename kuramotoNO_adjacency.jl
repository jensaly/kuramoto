using DifferentialEquations
using BenchmarkTools
using StaticArrays
using Plots
 
function kuramotoNO!(du, u, p, t)
    ω, K, N, uT, A1, A2, v1, v2 = p

    uT .= u'

    A1 .= uT .- u # Finding phase differences
    A2 .= sin.(A1)

    sum!(v1, A2) # Summing along sin(phase)
    mul!(v2, K, v1) # Adjancency matrix -> Interaction term

    du .= ω .+ v2 # Step

    nothing
end

N = 4

u = zeros(4)
du = similar(u)

ω = [6.6e9, 6.7e9, 6.2e9, 6.4e9]
K = fill(0.3e9, N, N)
K[diagind(K)] .= 0.0

u0 = zeros(N)
du = similar(u)

# Preallocate
uT = similar(u')
A1 = similar(u, N, N)
A2 = similar(u, N, N)
v1 = similar(u)
v2 = similar(u)

tstart = 0.0     # Start time
tend = 100e-9      # End time
dt = 1e-12        # Time step

tspan = (tstart, tend)  # Time span for simulation
p = (ω, K, N, uT, A1, A2, v1, v2)

# Define the ODE problem
prob = ODEProblem(kuramotoNO!, u0, tspan, p)

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