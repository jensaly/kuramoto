using DifferentialEquations
using BenchmarkTools
using StaticArrays
using Plots
using LinearAlgebra

mutable struct Kuramoto
    # Set elements
    u0::Vector{Float64}
    ω::Vector{Float64}
    K::Matrix{Float64}
    tstart::Float64
    tend::Float64
    dt::Float64
    
    # Private (preallocated)
    tspan::Tuple{Float64, Float64}
    du::Vector{Float64}
    A1::Matrix{Float64}
    A2::Matrix{Float64}
    v1::Vector{Float64}
    v2::Vector{Float64}

    Kuramoto(u0, ω, K, tstart, tend, dt) = new(u0, ω, K, tstart, tend, dt, (tstart, tend),
                                               zeros(N), zeros(N, N), zeros(N, N), zeros(N), zeros(N))
end

function run_kuramoto(model::Kuramoto)
    p = (model.ω, model.K, length(model.u0), (model.u0)', model.A1, model.A2, model.v1, model.v2)
    prob = ODEProblem(kuramotoNO!, model.u0, model.tspan, p)
    sol = solve(prob, Tsit5(), abstol=1e-10,reltol=1e-10, dt=model.dt)
end
 
function kuramotoNO!(du, u, p, t)
    ω, K, N, uT, A1, A2, v1, v2 = p

    A1 .= u' .- u # Finding phase differences
    A2 .= sin.(-A1)

    sum!(v1, A2) # Summing along sin(phase)
    mul!(v2, K, v1) # Adjancency matrix -> Interaction term

    du .= ω .+ v2 # Step

    nothing
end

function create_standard_K(k, N)
    K = fill(k, N, N) / N
    K[diagind(K)] .= 0.0
    return K
end

N = 4
K = 1e9

u = zeros(4)
du = similar(u)

ω = [6.6e9, 6.7e9, 6.2e9, 6.4e9]
K = create_standard_K(K, N)

tstart = 0.0     # Start time
tend = 100e-9      # End time
dt = 1e-12        # Time step

tspan = (tstart, tend)  # Time span for simulation
model = Kuramoto(u, ω, K, tstart, tend, dt)
@btime run_kuramoto(model)