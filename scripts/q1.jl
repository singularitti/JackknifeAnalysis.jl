using CSV
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

function hist_mean()
    for n in (1000, 10000)
        histograms = []
        for (i, variable) in enumerate(variables)
            X = sample(variable, PartitionSampler(n))
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
end

function stdmean()
    stds = hcat(map([1000, 10000]) do n
        map(variables) do variable
            round(truestd(variable, PartitionSampler(n)); digits=5)
        end
    end...)
    df = DataFrame(stds', [:v1, :v2, :v3, :v4, :v5])
    insertcols!(df, 1, :N => [1000, 10000])
    return latexify(df; env=:table)
end

function autocor_func()
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
    return savefig("tex/plots/Cn_C0.pdf")
end

function autocor_time()
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
    return latexify(df; env=:table)
end

minus(x, y) = x - y

function checkrelation(N=1000, nₘₐₓ=200)
    result = map(variables) do variable
        relation(variable, PartitionSampler(N), nₘₐₓ)
    end
    df = DataFrame(result)
    insertcols!(df, 1, :variable => [L"$v_%$i$" for i in 1:5])
    insertcols!(df, 1, :N => fill(N, 5))
    return combine(df, :, [:real, :estimated] => minus => :difference)
end

function checkrelations(N=[1000, 10000], nₘₐₓ=200)
    df = sort(vcat((checkrelation(n, nₘₐₓ) for n in N)...), :variable)
    CSV.write("relations.csv", df)
    return df
end

function covariance_matrix(variables)
    m = [cov(variable₁, variable₂) for variable₁ in variables, variable₂ in variables]
    return map(x -> round(x; digits=5), m)
end

function correlation_matrix(variables)
    ρ = [cor(variable₁, variable₂) for variable₁ in variables, variable₂ in variables]
    return map(x -> round(x; digits=5), ρ)
end

function sample_autocor_plot()
    plt = plot(;
        framestyle=:box,
        labelfontsize=12,
        tickfontsize=10,
        legendfontsize=12,
        right_margin=2mm,
    )
    for (i, variable) in enumerate(variables)
        s = sample(variable, PartitionSampler(10000))[1]
        C₀ = autocor(s, 0)
        N = 0:2000
        r = map(N) do n
            Cₙ = autocor(s, n)
            Cₙ / C₀
        end
        plot!(plt, N, r; label=L"$v_{%$i}$")
        xlims!(extrema(N))
        xlabel!(L"$n$")
        ylabel!(L"$C_{v_a}(n) / C_{v_a}(0)$")
    end
    return savefig("tex/plots/Cn_C0_sample.pdf")
end

function sample_autocor(N=[1000, 10000], nₘₐₓ=200)
    return map(N) do n
        map(variables) do variable
            s = sample(variable, PartitionSampler(n))[2]
            relation(s, nₘₐₓ)
        end
    end
end
