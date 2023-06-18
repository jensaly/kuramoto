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
    sol::ODESolution

    Kuramoto(u0, ω, K, tstart, tend, dt) = new(u0, ω, K, tstart, tend, dt, (tstart, tend),
                                               zeros(N), zeros(N, N), zeros(N, N), zeros(N), zeros(N))
end

function kuramotoNO_static!(du, u, p, t)
    ω, K, N, uT, A1, A2, v1, v2 = p

    A1 .= u' .- u # Finding phase differences
    A2 .= sin.(-A1)

    sum!(v1, A2) # Summing along sin(phase)
    mul!(v2, K, v1) # Adjancency matrix -> Interaction term

    du .= ω .+ v2 # Step

    nothing
end

function run_kuramoto_static(model::Kuramoto)
    p = (model.ω, model.K, length(model.u0), (model.u0)', model.A1, model.A2, model.v1, model.v2)
    prob = ODEProblem(kuramotoNO_static!, model.u0, model.tspan, p)
    model.sol = solve(prob, Tsit5(), abstol=1e-10,reltol=1e-10, dt=model.dt)
end

function kuramotoNO!(du, u, p, t)
    ω, K, N, uT, A1, A2, v1, v2 = p

    uT .= u'

    A1 .= uT .- u # Finding phase differences
    A2 .= sin.(-A1)

    sum!(v1, A2) # Summing along sin(phase)
    mul!(v2, K, v1) # Adjancency matrix -> Interaction term

    du .= ω .+ v2 # Step

    nothing
end

function run_kuramoto(model::Kuramoto)
    p = (model.ω, model.K, length(model.u0), (model.u0)', model.A1, model.A2, model.v1, model.v2)
    prob = ODEProblem(kuramotoNO!, model.u0, model.tspan, p)
    model.sol = solve(prob, Tsit5(), abstol=1e-10,reltol=1e-10, dt=model.dt)
end

function create_standard_K(k, N)
    K = fill(k, N, N) / N
    K[diagind(K)] .= 0.0
    return K
end

function plot_model_frequencies(model::Kuramoto)
    sol = model.sol
    t = sol.t
    θ1 = sol[1, :]
    θ2 = sol[2, :]
    θ3 = sol[3, :]
    θ4 = sol[4, :]

    f1 = diff(θ1) ./ diff(t) / 1e9
    f2 = diff(θ2) ./ diff(t) / 1e9
    f3 = diff(θ3) ./ diff(t) / 1e9
    f4 = diff(θ4) ./ diff(t) / 1e9

    p = plot(t[2:end], f1, label="Oscillator 1")
    plot!(p, t[2:end], f2, label="Oscillator 2")
    plot!(p, t[2:end], f3, label="Oscillator 3")
    plot!(p, t[2:end], f4, label="Oscillator 4")
    xlabel!(p, "Time")
    ylabel!(p, "Frequency")
    display(p)
end
