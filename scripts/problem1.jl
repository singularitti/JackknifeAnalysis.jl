using CSV: write
using DataFrames: DataFrame
using JackknifeAnalysis
using Latexify: latexify
using LaTeXStrings: @L_str
using Plots
using Plots.Measures

if !isdir("tex/plots")
    mkpath("tex/plots")
end

Plots.default(;
    framestyle=:box,
    labelfontsize=10,
    tickfontsize=8,
    legendfontsize=10,
    palette=:tab10,
    grid=nothing,
    legend_foreground_color=nothing,
)

const v1 = read("data/v1", Population)
const v2 = read("data/v2", Population)
const v3 = read("data/v3", Population)
const v4 = read("data/v4", Population)
const v5 = read("data/v5", Population)
const variables = [v1, v2, v3, v4, v5]

map(mean, (v1, v2, v3, v4, v5))

map(var, (v1, v2, v3, v4, v5))

function histogram_means()
    for n in (1000, 10000)
        histograms = []
        for (i, variable) in enumerate(variables)
            X = sampleby(variable, PartitionSampler(n))
            means = mean.(X)
            push!(
                histograms,
                histogram(
                    means;
                    legend=:none,
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

function plot_autocor(path="tex/plots/Cn_C0.pdf")
    plt = plot(; right_margin=2mm)
    for (i, variable) in enumerate(variables)
        C??? = autocor(variable, 0)
        N = 0:300
        r = map(N) do n
            C??? = autocor(variable, n)
            C??? / C???
        end
        plot!(plt, N, r; label=L"$v_{%$i}$")
        xlims!(extrema(N))
        xlabel!(L"$n$")
        ylabel!(L"$\hat{C}_{v_a}(n) / \hat{C}_{v_a}(0)$")
    end
    return savefig(path)
end

function plot_autocor_time(path="tex/plots/tau.pdf")
    plt = plot(; legend=:bottomright, right_margin=2mm)
    for (i, variable) in enumerate(variables)
        N = 1:200
        ????? = map(Base.Fix1(int_autocor_time, variable), N)
        plot!(plt, N, ?????; label=L"$v_{%$i}$")
        xlims!(extrema(N))
        xlabel!(L"$n_\textnormal{cut}$")
        ylabel!(L"$\hat{\tau}_{v_a}$")
    end
    savefig(path)
    ????? = map(variables) do variable
        round(int_autocor_time(variable, 200); digits=5)
    end
    df = DataFrame(?????', [:v1, :v2, :v3, :v4, :v5])
    return latexify(df; env=:table)
end

minus(x, y) = x - y

function checkrelation(N=1000, n?????????=200)
    result = map(variables) do variable
        relation(variable, PartitionSampler(N), n?????????)
    end
    df = DataFrame(result)
    insertcols!(df, 1, :variable => [L"$v_%$i$" for i in 1:5])
    insertcols!(df, 1, :N => fill(N, 5))
    return combine(df, :, [:real, :estimated] => minus => :difference)
end

function checkrelations(N=[1000, 10000], n?????????=200)
    df = sort(vcat((checkrelation(n, n?????????) for n in N)...), :variable)
    write("relations.csv", df)
    return df
end

function covariance_matrix(variables)
    m = [cov(variable???, variable???) for variable??? in variables, variable??? in variables]
    return map(x -> round(x; digits=5), m)
end

function correlation_matrix(variables)
    ?? = [cor(variable???, variable???) for variable??? in variables, variable??? in variables]
    return map(x -> round(x; digits=5), ??)
end

function plot_sample_autocor(N, n?????????=20, index=1, path="tex/plots/Cn_C0_sample_$N.pdf")
    plot(; right_margin=2mm)
    hline!([0]; color="black", label="")
    for (i, variable) in enumerate(variables)
        s = sampleby(variable, PartitionSampler(N))[index]
        C??? = autocor(s, 0)
        terms = 0:n?????????
        r = map(terms) do n
            C??? = autocor(s, n)
            C??? / C???
        end
        plot!(terms, r; label=L"$v_{%$i}$")
        xlims!(extrema(terms))
        xlabel!(L"$n$")
        ylabel!(L"$C_{v_a}(n) / C_{v_a}(0)$")
    end
    return savefig(path)
end

function plot_autocor_time(N, n?????????=200, index=1, path="tex/plots/tau_sample_$N.pdf")
    plt = plot(; legend=:bottom, right_margin=2mm)
    for (i, variable) in enumerate(variables)
        sample = sampleby(variable, PartitionSampler(N))[index]
        ????? = map(Base.Fix1(int_autocor_time, sample), 1:n?????????)
        plot!(plt, 1:n?????????, ?????; label=L"$v_{%$i}$")
        xlims!(extrema(1:n?????????))
        xlabel!(L"$n_\textnormal{cut}$")
        ylabel!(L"$\tau_{v_a}$")
    end
    return savefig(path)
end

function guess_n?????????(variable, N, n?????????=200, index=1)
    s = sampleby(variable, PartitionSampler(N))[index]
    terms = 1:n?????????
    for n in terms
        C??? = autocor(s, n)
        if C??? < 0
            return n - 1
        end
    end
end

getsample(variable, N=1000, index=1) = sampleby(variable, PartitionSampler(N))[index]

function get_sample_relation(N, index=1)
    return map(variables) do variable
        sample = getsample(variable, N, index)
        n????????? = guess_n?????????(variable, N, 400, index)
        relation(sample, n?????????)
    end
end

function sample_correlation_matrix(samples)
    ?? = [cor(sample???, sample???) for sample??? in samples, sample??? in samples]
    return map(x -> round(x; digits=5), ??)
end
