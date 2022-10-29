using JackknifeAnalysis
using LaTeXStrings
using Plots
using Plots.Measures

if !isdir("tex/plots")
    mkpath("tex/plots")
end

Plots.default(;
    framestyle=:box,
    labelfontsize=12,
    tickfontsize=10,
    legendfontsize=12,
    palette=:tab20,
    grid=nothing,
    legend_foreground_color=nothing,
)

const v₁ = read("data/v1", Population)
const v₂ = read("data/v2", Population)
const v₃ = read("data/v3", Population)
const v₄ = read("data/v4", Population)
const v₅ = read("data/v5", Population)
const variables = [v₁, v₂, v₃, v₄, v₅]
const 𝐛 = [2, 4, 5, 8, 10, 20, 25, 40, 50, 100, 125, 200, 250]

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

function getsamples(variable::Population, samplesize=5000)
    return sampleby(variable, PartitionSampler(samplesize))
end

function sample_means(variable::Population, samplesize=5000)
    return mean.(sampleby(variable, PartitionSampler(samplesize)))
end

function sample_truemean(variable::Population, samplesize=5000)
    return mean(sample_means(variable, samplesize))
end

function sample_variance(variable::Population, samplesize=5000)
    return var(variable, PartitionSampler(samplesize))
end

function var_f₁(samplesize=5000)
    𝕍₁, 𝕍₂ = sample_variance(v₁, samplesize), sample_variance(v₂, samplesize)
    v̄̂₁, v̄̂₂ = sample_truemean(v₁, samplesize), sample_truemean(v₂, samplesize)
    return 1 / v̄̂₂^2 * 𝕍₁ + v̄̂₁^2 / v̄̂₂^4 * 𝕍₂
end

function var_f₂(samplesize=5000)
    𝕍₃, 𝕍₄ = sample_variance(v₃, samplesize), sample_variance(v₄, samplesize)
    v̄̂₃, v̄̂₄ = sample_truemean(v₃, samplesize), sample_truemean(v₄, samplesize)
    return exp(v̄̂₃ - v̄̂₄)^2 * (𝕍₃ + 𝕍₄)
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

function plot_f₁_f₂(samplesize=5000)
    plot(; xlabel="samples", ylabel=L"$f_i(\bar{v}_a)$", legend=:right)
    results = f₁(v₁, v₂, samplesize)
    scatter!(
        eachindex(results),
        results;
        label=L"$f_1(\bar{v}_1, \bar{v}_2)$",
        markersizes=2,
        markerstrokewidth=0,
    )
    plot!(
        eachindex(results),
        results;
        label="",
        xlims=(firstindex(results), lastindex(results)),
    )
    results = f₂(v₃, v₄, samplesize)
    scatter!(
        eachindex(results),
        results;
        label=L"$f_2(\bar{v}_3, \bar{v}_4)$",
        markersizes=2,
        markerstrokewidth=0,
    )
    plot!(
        eachindex(results),
        results;
        label="",
        xlims=(firstindex(results), lastindex(results)),
    )
    return savefig("tex/plots/f1_f2.pdf")
end

truemean_f₁(samplesize=5000) = mean(f₁(v₁, v₂, samplesize))

truestd_f₁(samplesize=5000) = std(Sample(f₁(v₁, v₂, samplesize)))

truemean_f₂(samplesize=5000) = mean(f₂(v₃, v₄, samplesize))

truestd_f₂(samplesize=5000) = std(Sample(f₂(v₃, v₄, samplesize)))

function guess_nₘₐₓ(variable, nₘₐₓ=200, i=1)
    s = sampleby(variable, PartitionSampler(5000))[i]
    terms = 1:nₘₐₓ
    for n in terms
        Cₙ = autocor(s, n)
        if Cₙ < 0
            return n - 1
        end
    end
end

function getbins(variable, i, samplesize=5000, b=200)
    sample = sampleby(variable, PartitionSampler(samplesize))[i]
    return Sample(mean.(sampleby(Population(sample), PartitionSampler(b))))
end

jackknife(sample) = sampleby(sample, JackknifeSampler())

function plot_jackknife_means(i=1)
    plot(;
        xlims=extrema(𝐛),
        xlabel=L"size of bin ($b$)",
        ylabel=L"$\bar{v}_a$",
        legend_column=2,
    )
    for (j, variable) in enumerate(variables)
        means = map(𝐛) do b
            mean(jackknife(getbins(variable, i, 5000, b)))
        end
        scatter!(𝐛, means; label=L"$\bar{v}_{%$j}$", markersizes=2, markerstrokewidth=0)
        plot!(𝐛, means; label="")
    end
    savefig("tex/plots/JA_sample_means.pdf")
    plot(;
        xlims=extrema(𝐛),
        yformatter=x -> "$(round(x/1e-15; digits=3))",
        xlabel=L"size of bin ($b$)",
        ylabel=L"$(\bar{v}_a - \hat{\bar{v}}_a) / 10^{-15}$",
    )
    for (j, variable) in enumerate(variables)
        means = map(𝐛) do b
            mean(jackknife(getbins(variable, i, 5000, b)))
        end
        μ = mean(means)
        deviations = map(means) do mean
            mean - μ
        end
        scatter!(
            𝐛,
            deviations;
            label=L"$\bar{v}_{%$i} - \hat{\bar{v}}_{%$j}$",
            markersizes=2,
            markerstrokewidth=0,
        )
        plot!(𝐛, deviations; label="")
    end
    savefig("tex/plots/JA_sample_deviations.pdf")
    return nothing
end

function plot_jackknife_std(i=1)
    plot(;
        xlims=extrema(𝐛),
        xlabel=L"size of bin ($b$)",
        ylabel=L"$\sigma_{\bar{v}_a}$",
        legend=:bottomleft,
    )
    for (j, variable) in enumerate(variables)
        stds = map(𝐛) do b
            std(jackknife(getbins(variable, i, 5000, b)))
        end
        scatter!(
            𝐛,
            stds;
            ylims=(0, Inf),
            label=L"$\sigma_{\bar{v}_{%$j}}$",
            markersizes=2,
            markerstrokewidth=0,
        )
        plot!(𝐛, stds; label="")
    end
    savefig("tex/plots/JA_sample_std.pdf")
    return nothing
end

function getsample(variable₁, variable₂, i, b)
    return map(variable -> jackknife(getbins(variable, i, 5000, b)), (variable₁, variable₂))
end

function get_f_mean(f, variable₁, variable₂, i=1)
    return map(𝐛) do b
        sample₁, sample₂ = getsample(variable₁, variable₂, i, b)
        f(mean(sample₁), mean(sample₂))
    end
end

function get_jackknife_f_mean(f, variable₁, variable₂, i, b)
    sample₁, sample₂ = getsample(variable₁, variable₂, i, b)
    return mean([f(x′, y′) for (x′, y′) in zip(sample₁, sample₂)])
end
function get_jackknife_f_mean(::typeof(f₁), i=1)
    return map(𝐛) do b
        get_jackknife_f_mean(f₁, v₁, v₂, i, b)
    end
end
function get_jackknife_f_mean(::typeof(f₂), i=1)
    return map(𝐛) do b
        get_jackknife_f_mean(f₂, v₃, v₄, i, b)
    end
end

function plot_jackknife_f_mean(i=1)
    plt = plot(;
        xlims=extrema(𝐛),
        xlabel=L"size of bin ($b$)",
        ylabel=L"$\Delta \bar{f}_i / 10^{-5}$",
        legend=:right,
    )
    means = get_jackknife_f_mean(f₁, i)
    f̄ = get_f_mean(f₁, v₁, v₂, i)
    data = (means .- f̄) / 1e-5
    scatter!(plt, 𝐛, data; label=L"$f_1$", markersizes=2, markerstrokewidth=0)
    plot!(plt, 𝐛, data; label="")
    means = get_jackknife_f_mean(f₂, i)
    f̄ = get_f_mean(f₂, v₃, v₄, i)
    data = (means .- f̄) / 1e-5
    scatter!(plt, 𝐛, data; label=L"$f_2$", markersizes=2, markerstrokewidth=0)
    plot!(plt, 𝐛, data; label="")
    savefig("tex/plots/JA_f1_f2_mean.pdf")
    return nothing
end

function get_f_std(f, variable₁, variable₂, i, b)
    sample₁, sample₂ = map(
        variable -> jackknife(getbins(variable, i, 5000, b)), (variable₁, variable₂)
    )
    return std(f, sample₁, sample₂)
end
function get_f_std(::typeof(f₁), i=1)
    return map(𝐛) do b
        get_f_std(f₁, v₁, v₂, i, b)
    end
end
function get_f_std(::typeof(f₂), i=1)
    return map(𝐛) do b
        get_f_std(f₂, v₃, v₄, i, b)
    end
end

function plot_jackknife_f_std(::typeof(f₁), i=1)
    plot(;
        xlims=extrema(𝐛),
        xlabel=L"size of bin ($b$)",
        ylabel=L"$\sigma_{f_1,N}$",
        legend=:none,
    )
    stds = get_f_std(f₁, i)
    scatter!(𝐛, stds; ylims=extrema(stds), markersizes=2, markerstrokewidth=0)
    plot!(𝐛, stds; label="")
    savefig("tex/plots/JA_f1_std.pdf")
    return nothing
end
function plot_jackknife_f_std(::typeof(f₂), i=1)
    plot(;
        xlims=extrema(𝐛),
        xlabel=L"size of bin ($b$)",
        ylabel=L"$\sigma_{f_2,N}$",
        legend=:none,
    )
    stds = get_f_std(f₂, i)
    scatter!(𝐛, stds; ylims=extrema(stds), markersizes=2, markerstrokewidth=0)
    plot!(𝐛, stds; label="")
    savefig("tex/plots/JA_f2_std.pdf")
    return nothing
end
