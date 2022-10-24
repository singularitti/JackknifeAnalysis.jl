using JackknifeAnalysis
using LaTeXStrings
using Plots
using Plots.Measures

function initplot()
    if !isdir("tex/plots")
        mkpath("tex/plots")
    end
    return Plots.default(;
        framestyle=:box,
        labelfontsize=12,
        tickfontsize=10,
        legendfontsize=12,
        palette=:tab20,
        grid=nothing,
        legend_foreground_color=nothing,
    )
end

initplot()

const v1 = read("data/v1", Population)
const v2 = read("data/v2", Population)
const v3 = read("data/v3", Population)
const v4 = read("data/v4", Population)
const v5 = read("data/v5", Population)
const variables = [v1, v2, v3, v4, v5]

f₁(v̄₁, v̄₂) = v̄₁ / v̄₂
f₁(sample₁::Sample, sample₂::Sample) = f₁(mean(sample₁), mean(sample₂))
function f₁(v1::Population, v2::Population, samplesize=5000)
    sampler = PartitionSampler(samplesize)
    return map(sampleby(v1, sampler), sampleby(v2, sampler)) do sample₁, sample₂
        f₁(sample₁, sample₂)
    end
end

f₂(v̄₃, v̄₄) = exp(v̄₃ - v̄₄)
f₂(sample₃::Sample, sample₄::Sample) = f₂(mean(sample₃), mean(sample₄))
function f₂(v3::Population, v4::Population, samplesize=5000)
    sampler = PartitionSampler(samplesize)
    return map(sampleby(v3, sampler), sampleby(v4, sampler)) do sample₃, sample₄
        f₂(sample₃, sample₄)
    end
end

function sample_means(variable::Population, samplesize=5000)
    X = sampleby(variable, PartitionSampler(samplesize))
    return mean.(X)
end

function plot_sample_means(samplesize=5000)
    plot(; xlabel="samples", ylabel=L"$\bar{v}_a$")
    for (i, variable) in enumerate(variables)
        means = sample_means(variable, samplesize)
        scatter!(
            eachindex(means),
            means;
            label=L"$\bar{v}_{%$i, N}$",
            markersizes=2,
            markerstrokewidth=0,
            legend_columns=5,
        )
        plot!(
            eachindex(means), means; label="", xlims=(firstindex(means), lastindex(means))
        )
    end
    return savefig("tex/plots/sample_means.pdf")
end

function plot_f₁(samplesize=5000)
    plot(; legend=:none, xlabel="samples", ylabel=L"$f_1(\bar{v}_1, \bar{v}_2)$")
    results = f₁(v1, v2, samplesize)
    scatter!(eachindex(results), results; markersizes=2, markerstrokewidth=0)
    plot!(eachindex(results), results; xlims=(firstindex(results), lastindex(results)))
    return savefig("tex/plots/f1.pdf")
end

function plot_f₂(samplesize=5000)
    plot(; legend=:none, xlabel="samples", ylabel=L"$f_2(\bar{v}_3, \bar{v}_4)$")
    results = f₂(v3, v4, samplesize)
    scatter!(eachindex(results), results; markersizes=2, markerstrokewidth=0)
    plot!(eachindex(results), results; xlims=(firstindex(results), lastindex(results)))
    return savefig("tex/plots/f2.pdf")
end

truemean_f₁(samplesize=5000) = mean(f₁(v1, v2, samplesize))

truestd_f₁(samplesize=5000) = std(Sample(f₁(v1, v2, samplesize)))

truemean_f₂(samplesize=5000) = mean(f₂(v3, v4, samplesize))

truestd_f₂(samplesize=5000) = std(Sample(f₂(v3, v4, samplesize)))

function getbins(variable, index, samplesize=5000, binsize=200)
    sample = sampleby(variable, PartitionSampler(samplesize))[index]
    return Sample(mean.(sampleby(Population(sample), PartitionSampler(binsize))))
end

jackknife(sample) = sampleby(sample, JackknifeSampler())

function plot_jackknife_means(index=1)
    binsizes = [
        2, 4, 5, 8, 10, 20, 25, 40, 50, 100, 125, 200, 250, 500, 625, 1000, 1250, 2500
    ]
    plot(;
        xlims=extrema(binsizes),
        xlabel=L"size of bin ($b$)",
        ylabel=L"$\bar{v}_a$",
        legend_column=2,
    )
    for (i, variable) in enumerate(variables)
        means = map(binsizes) do binsize
            mean(jackknife(getbins(variable, index, 5000, binsize)))
        end
        scatter!(
            binsizes, means; label=L"$\bar{v}_{%$i}$", markersizes=2, markerstrokewidth=0
        )
        plot!(binsizes, means; label="")
    end
    savefig("tex/plots/JA_sample_means.pdf")
    plot(;
        xlims=extrema(binsizes),
        yformatter=x -> "$(round(x/1e-15; digits=3))",
        xlabel=L"size of bin ($b$)",
        ylabel=L"$(\bar{v}_a - \hat{\bar{v}}_a) / 10^{-15}$",
    )
    for (i, variable) in enumerate(variables)
        means = map(binsizes) do binsize
            mean(jackknife(getbins(variable, index, 5000, binsize)))
        end
        μ = mean(means)
        deviations = map(means) do mean
            mean - μ
        end
        scatter!(
            binsizes,
            deviations;
            label=L"$\bar{v}_{%$i} - \hat{\bar{v}}_{%$i}$",
            markersizes=2,
            markerstrokewidth=0,
        )
        plot!(binsizes, deviations; label="")
    end
    savefig("tex/plots/JA_sample_deviations.pdf")
    return nothing
end

function plot_jackknife_std(index=1)
    binsizes = [
        2, 4, 5, 8, 10, 20, 25, 40, 50, 100, 125, 200, 250, 500, 625, 1000, 1250, 2500
    ]
    plot(;
        xlims=extrema(binsizes),
        xlabel=L"size of bin ($b$)",
        ylabel=L"$\sigma_{\bar{v}_a}$",
        legend=:topleft,
    )
    for (i, variable) in enumerate(variables)
        stds = map(binsizes) do binsize
            std(jackknife(getbins(variable, index, 5000, binsize)))
        end
        scatter!(
            binsizes,
            stds;
            ylims=(0, Inf),
            label=L"$\sigma_{\bar{v}_{%$i}}$",
            markersizes=2,
            markerstrokewidth=0,
        )
        plot!(binsizes, stds; label="")
    end
    savefig("tex/plots/JA_sample_std.pdf")
    return nothing
end
