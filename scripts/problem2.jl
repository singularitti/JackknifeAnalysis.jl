using JackknifeAnalysis
using LaTeXStrings
using Plots
using Plots.Measures

if !isdir("tex/plots")
    mkpath("tex/plots")
end

v1 = read("data/v1", Population)
v2 = read("data/v2", Population)
v3 = read("data/v3", Population)
v4 = read("data/v4", Population)
v5 = read("data/v5", Population)
variables = [v1, v2, v3, v4, v5]

f₁(v̄₁, v̄₂) = v̄₁ / v̄₂
f₁(sample₁::Sample, sample₂::Sample) = f₁(mean(sample₁), mean(sample₂))
function f₁(v1::Population, v2::Population, samplesize=5000)
    sampler = PartitionSampler(samplesize)
    return map(sample(v1, sampler), sample(v2, sampler)) do sample₁, sample₂
        f₁(sample₁, sample₂)
    end
end

f₂(v̄₃, v̄₄) = exp(v̄₃ - v̄₄)
f₂(sample₃::Sample, sample₄::Sample) = f₂(mean(sample₃), mean(sample₄))
function f₂(v3::Population, v4::Population, samplesize=5000)
    sampler = PartitionSampler(samplesize)
    return map(sample(v3, sampler), sample(v4, sampler)) do sample₃, sample₄
        f₂(sample₃, sample₄)
    end
end

function sample_means(variable::Population, samplesize=5000)
    X = sampleby(variable, PartitionSampler(samplesize))
    return mean.(X)
end

function plot_sample_means(samplesize=5000)
    plot(;
        framestyle=:box,
        xlabel="samples",
        ylabel=L"$\bar{v}_a$",
        legend=:bottom,
        labelfontsize=12,
        tickfontsize=10,
        legendfontsize=12,
        palette=:tab20,
    )
    for (i, variable) in enumerate(variables)
        means = sample_means(variable, samplesize)
        scatter!(eachindex(means), means; label="", markersizes=2, markerstrokewidth=0)
        plot!(
            eachindex(means),
            means;
            label=L"$\bar{v}_{%$i, N}$",
            xlims=(firstindex(means), lastindex(means)),
        )
    end
    return savefig("tex/plots/JA_sample_means.pdf")
end

function plot_f₁(samplesize=5000)
    plot(;
        framestyle=:box,
        legend=:none,
        xlabel="samples",
        ylabel=L"$f_1(\bar{v}_1, \bar{v}_2)$",
        tickfontsize=10,
        palette=:tab20,
    )
    results = f₁(v1, v2, samplesize)
    scatter!(eachindex(results), results; markersizes=2, markerstrokewidth=0)
    plot!(eachindex(results), results; xlims=(firstindex(results), lastindex(results)))
    return savefig("tex/plots/f1.pdf")
end

function plot_f₂(samplesize=5000)
    plot(;
        framestyle=:box,
        legend=:none,
        xlabel="samples",
        ylabel=L"$f_2(\bar{v}_3, \bar{v}_4)$",
        tickfontsize=10,
        palette=:tab20,
    )
    results = f₂(v3, v4, samplesize)
    scatter!(eachindex(results), results; markersizes=2, markerstrokewidth=0)
    plot!(eachindex(results), results; xlims=(firstindex(results), lastindex(results)))
    return savefig("tex/plots/f2.pdf")
end

truemean_f₁(samplesize=5000) = mean(f₁(v1, v2, samplesize))

truestd_f₁(samplesize=5000) = std(Sample(f₁(v1, v2, samplesize)))

truemean_f₂(samplesize=5000) = mean(f₂(v3, v4, samplesize))

truestd_f₂(samplesize=5000) = std(Sample(f₂(v3, v4, samplesize)))
