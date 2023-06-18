include("kuramotoNO_dynamic_struct.jl")

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
K = 0.97e9

u = zeros(4)

ω = [6.6e9, 6.7e9, 6.2e9, 6.4e9]
K = create_standard_K(K, N) * N

K[1,2] /= 3
K[2,1] /= 3
K[3,4] /= 3
K[4,3] /= 3

tstart = 0.0     # Start time
tend = 100e-9      # End time
dt = 1e-12        # Time step

tspan = (tstart, tend)  # Time span for simulation
model = Kuramoto(u, ω, K, tstart, tend, dt)

#display(K)

freq_diff = zeros(300,300)
natfreq1 = range(96e9, 104e9, 300)
natfreq2 = range(96e9, 104e9, 300)

#display(natfreq1)


for i in 1:length(natfreq1)
    for j in 1:length(natfreq2)
        model.ω = [natfreq1[i], natfreq2[j], 99.5e9, 100.5e9]
        run_kuramoto(model)
        freq_diff[i,j] = output_freq_diff(model)
        #print(model.ω, "\n")
    end
    print(i)
end

filename = "matrix_output.txt"

writedlm(filename, freq_diff, '\t')
print("Exit code 0")
#=
run_kuramoto(model)

plot_model_frequencies(model)

model.ω = [6.6e9, 6.7e9, 6.2e9, 4.4e9]

run_kuramoto(model)

plot_model_frequencies(model)

model.ω = [6.6e9, 6.7e9, 6.2e9, 2.4e9]

run_kuramoto(model)

plot_model_frequencies(model)
=#