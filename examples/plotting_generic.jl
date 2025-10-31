using Plots

function plot_model_frequencies(model::KuramotoModel, step::Int64=1, legend::Bool=true)
    sol = model.sol
    t = sol.t
    p = plot(legend=false)
    if legend
        for i in 1:length(model.ω)
            plot!(p, t[2:step:end], (diff(sol[i, :]) ./ diff(t))[1:step:end] / 1e9, label="Oscillator " * string(i), linewidth=5)
        end
    else
        for i in 1:length(model.ω)
            plot!(p, t[2:step:end], (diff(sol[i, :]) ./ diff(t))[1:step:end] / 1e9, linewidth=5)
        end
    end
    xlabel!(p, "Time")
    ylabel!(p, "Frequency")
end

function plot_model_phases(model::KuramotoModel)
    sol = model.sol
    t = sol.t
    p = plot()
    for i in 1:length(model.ω)
        plot!(p, t, sol[i, :], label="Oscillator " * string(i))
    end
    xlabel!(p, "Time")
    ylabel!(p, "Frequency")
    #display(p)
end

function order_parameter(model::KuramotoModel)
    sol = model.sol
    N = length(model.u0)
    cosmat = cos.(model.sol)
    sinmat = sin.(model.sol)
    sum_ = similar(cosmat)
    fill!(sum_, 0.0)
    for i in 1:N
        for j in 1:N
            if i != j
                sum_[i,:] .= cosmat[i,:].*cosmat[j,:] .+ sinmat[i,:].*sinmat[j,:]
            end
        end
    end
    sum_ = sqrt.(N .+ sum(sum_, dims=1)) ./ N
    return sum_
end