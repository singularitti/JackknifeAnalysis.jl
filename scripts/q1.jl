using DataFrames
using JackknifeAnalysis
using Latexify
using LaTeXStrings
using Plots
using Plots.Measures

if !isdir("tex/plots")
    mkpath("tex/plots")
end

const v1 = read("data/v1", Population)
const v2 = read("data/v2", Population)
const v3 = read("data/v3", Population)
const v4 = read("data/v4", Population)
const v5 = read("data/v5", Population)
const variables = [v1, v2, v3, v4, v5]

map(mean, (v1, v2, v3, v4, v5))

map(var, (v1, v2, v3, v4, v5))

let
    for n in (1000, 10000)
        histograms = []
        for (i, variable) in enumerate(variables)
            X = sample(variable, ContinuousSampler(n))
            means = mean.(X)
            push!(
                histograms,
                histogram(
                    means;
                    legend=:none,
                    framestyle=:box,
                    xlims=extrema(means),
                    ylims=(0, Inf),
                    xlabel=L"$\bar{v}_%$i$",
                    ylabel=i in (1, 4) ? "frequency" : "",
                    labelfontsize=8,
                    tickfontsize=5,
                    top_margin=-1mm,
                    bottom_margin=-1mm,
                    left_margin=0mm,
                    right_margin=-1mm,
                ),
            )
            plot(histograms...)
        end
        savefig("tex/plots/sample means N=$n.pdf")
    end
    stds = hcat(map([1000, 10000]) do n
        map(variables) do variable
            round(truestd(variable, ContinuousSampler(n)); digits=5)
        end
    end...)
    df = DataFrame(stds', [:v1, :v2, :v3, :v4, :v5])
    insertcols!(df, 1, :N => [1000, 10000])
    latexify(df; env=:table)
end

let
    plt = plot(;
        framestyle=:box,
        labelfontsize=12,
        tickfontsize=10,
        legendfontsize=12,
        right_margin=2mm,
    )
    for (i, variable) in enumerate(variables)
        C₀ = truecor(variable, 0)
        N = 0:300
        r = map(N) do n
            Cₙ = truecor(variable, n)
            Cₙ / C₀
        end
        plot!(plt, N, r; label=L"$v_{%$i}$")
        xlims!(extrema(N))
        xlabel!(L"$n$")
        ylabel!(L"$\hat{C}_{v_a}(n) / \hat{C}_{v_a}(0)$")
    end
    savefig("tex/plots/Cn_C0.pdf")
end

let
    plt = plot(;
        framestyle=:box,
        legend=:bottomright,
        labelfontsize=12,
        tickfontsize=10,
        legendfontsize=12,
        right_margin=2mm,
    )
    for (i, variable) in enumerate(variables)
        N = 1:200
        τₙ = map(Base.Fix1(int_autocor_time, variable), N)
        plot!(plt, N, τₙ; label=L"$v_{%$i}$")
        xlims!(extrema(N))
        xlabel!(L"$n_\textnormal{cut}$")
        ylabel!(L"$\hat{\tau}_{v_a}$")
    end
    savefig("tex/plots/tau.pdf")
    τᵥ = map(variables) do variable
        round(int_autocor_time(variable, 200); digits=5)
    end
    df = DataFrame(τᵥ', [:v1, :v2, :v3, :v4, :v5])
    latexify(df; env=:table)
end
