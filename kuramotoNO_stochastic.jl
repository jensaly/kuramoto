using DifferentialEquations
using BenchmarkTools
using StaticArrays
using Plots
using LinearAlgebra

 
function kuramotoNO(u, p, t)
    ω, K, N, uT, A1, A2, v1, v2, nV = p

    A1 = u' .- u # Phase differences
    A2 = sin.(-A1)

    v1 = sum(A2, dims=2) # Summing along rows
    v2 = K * v1 # Adjancency matrix -> Interaction term

    return SVector{4, Float64}(ω .+ v2)
end

function kuramotoNO_noise(u, p, t)
    ω, K, N, uT, A1, A2, v1, v2, nV = p
    return nV # Uncorrelated white noise of power D
end

N = 4
K = 0.5e9
D = 0.000005e9
W = WienerProcess(0.0,0.0,0.0)

u = @SVector zeros(N)
du = @SVector zeros(N)

ω = @SVector [6.6e9, 6.7e9, 6.2e9, 6.4e9]
K = fill(K, N, N) / N
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
nV = @SVector fill(D, N)
nV = nV ./ N # Normalizing the noise

# ODE setup
tstart = 0.0     # Start time
tend = 100e-9      # End time
dt = 1e-15       # Time step

tspan = (tstart, tend)  # Time span for simulation
p = (ω, K, N, uT, A1, A2, v1, v2, nV)

# Define the ODE problem
prob = SDEProblem(kuramotoNO, kuramotoNO_noise, u0, tspan, p, noise=W)

# Solve the ODE problem
#@btime sol = solve(prob, Tsit5(), abstol=1e-10,reltol=1e-10, dt=dt)
sol = solve(prob, abstol=1e-7,reltol=1e-7, dt=dt)


# Access the solution
t = sol.t
θ1 = sol[1, :]
θ2 = sol[2, :]
θ3 = sol[3, :]
θ4 = sol[4, :]

display(t)

f1 = diff(θ1) ./ diff(t) / 1e9
f2 = diff(θ2) ./ diff(t) / 1e9
f3 = diff(θ3) ./ diff(t) / 1e9
f4 = diff(θ4) ./ diff(t) / 1e9

ind = trunc(Int, length(t)/2)+1
print(length(t))

plot(t[ind+1:end], f1[ind:end], label="Oscillator 1")
plot!(t[ind+1:end], f2[ind:end], label="Oscillator 2")
plot!(t[ind+1:end], f3[ind:end], label="Oscillator 3")
plot!(t[ind+1:end], f4[ind:end], label="Oscillator 4")
xlabel!("Time")
ylabel!("Frequency")
