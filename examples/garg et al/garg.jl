include("../../kuramoto.jl")

using Distributed
using DelimitedFiles

function output_freq_diff(model::Kuramoto)
    sol = model.sol
    t = sol.t
    θ1 = sol[1, :]
    θ2 = sol[2, :]
    θ3 = sol[3, :]
    θ4 = sol[4, :]

    idx = findfirst(x -> x > 70e-9, t)  # Find the first index where t > threshold

    f1 = diff(θ1) ./ diff(t) / 1e9
    f2 = diff(θ2) ./ diff(t) / 1e9
    f3 = diff(θ3) ./ diff(t) / 1e9
    f4 = diff(θ4) ./ diff(t) / 1e9

    return abs(mean(f3[idx:end]) - mean(f4[idx:end]))
end

N = 4
K = 0.05e9

u = zeros(4)

ω = [0,0,0,0] # Placeholder
K = create_standard_K(K, N)

K[1,2] = 0.0398e9 / 4
K[2,1] = 0.0398e9 / 4
K[3,4] = 0.0398e9 / 4
K[4,3] = 0.0398e9 / 4

K = 4 * K
display(K)

tstart = 0.0     # Start time
tend = 100e-9      # End time
dt = 1e-12        # Time step

tspan = (tstart, tend)  # Time span for simulation
model = Kuramoto(u, ω, K, tstart, tend, dt)

freq_diff = zeros(300,300)
natfreq1 = range(6.20e9, 7.20e9, 300)
natfreq2 = range(6.20e9, 7.20e9, 300)

for i in 1:length(natfreq1)
    for j in 1:length(natfreq2)
        model.ω = [natfreq1[i], natfreq2[j], 6.63e9, 6.81e9]
        run_kuramoto(model, 1e-12, 1e-12)
        freq_diff[i,j] = output_freq_diff(model)
    end
    print(i)
end

filename = "examples/garg et al/garg.txt"

writedlm(filename, freq_diff, '\t')
print("Exit code 0")

