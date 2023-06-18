using DifferentialEquations
using BenchmarkTools
using StaticArrays
using Plots
using LinearAlgebra
 
function kuramoto!(du, u, p, t)
    ω, K, N, uT, A1, A2, v1, v2 = p

    uT .= u'

    A1 .= uT .- u # Finding phase differences
    A2 .= sin.(-A1)

    sum!(v1, A2) # Summing along sin(phase)
    mul!(v2, K, v1) # Adjancency matrix -> Interaction term

    du .= ω .+ v2 # Step

    nothing
end

N = 4
K = 1e9

u = zeros(4)
du = similar(u)

ω = [6.6e9, 6.7e9, 6.2e9, 6.4e9]
K = fill(K, N, N) / N
K[diagind(K)] .= 0.0
display(typeof(K))

u0 = zeros(N)
du = similar(u)

# Preallocate
uT = similar(u')
A1 = similar(u, N, N)
A2 = similar(u, N, N)
v1 = similar(u)
v2 = similar(u)

display(u)
display(uT)

tstart = 0.0     # Start time
tend = 100e-9      # End time
dt = 1e-12        # Time step

tspan = (tstart, tend)  # Time span for simulation
p = (ω, K, N, uT, A1, A2, v1, v2)

# Define the ODE problem
prob = ODEProblem(kuramoto!, u0, tspan, p)

# Solve the ODE problem
@btime sol = solve(prob, Tsit5(), abstol=1e-10,reltol=1e-10, dt=dt)

#=
# Access the solution
t = sol.t
θ1 = sol[1, :]
θ2 = sol[2, :]
θ3 = sol[3, :]
θ4 = sol[4, :]

f1 = diff(θ1) ./ diff(t) / 1e9
f2 = diff(θ2) ./ diff(t) / 1e9
f3 = diff(θ3) ./ diff(t) / 1e9
f4 = diff(θ4) ./ diff(t) / 1e9

plot(t[2:end], f1, label="Oscillator 1")
plot!(t[2:end], f2, label="Oscillator 2")
plot!(t[2:end], f3, label="Oscillator 3")
plot!(t[2:end], f4, label="Oscillator 4")
xlabel!("Time")
ylabel!("Frequency")
=#