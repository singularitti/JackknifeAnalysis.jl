using JackknifeAnalysis
using LaTeXStrings: @L_str
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
    palette=:tab10,
    grid=nothing,
    legend_foreground_color=nothing,
)

function plot_autocor(
    u::Sample, k::Sample, nₘₐₓ=3000, path="tex/plots/MD_autocor_energy.pdf"
)
    plot(; right_margin=3mm)
    xlabel!(L"$n$")
    ylabel!(L"$\hat{C}(n) / \hat{C}(0)$")
    for (sample, label) in
        zip((u, k, u .+ k), ("potential", "kinetic", "total") .* " energy")
        C₀ = autocor(sample, 0)
        N = 0:nₘₐₓ
        r = map(N) do n
            Cₙ = autocor(sample, n)
            Cₙ / C₀
        end
        plot!(N, r; label=label)
        xlims!(extrema(N))
    end
    return savefig(path)
end

function guess_nₘₐₓ(sample)
    return findfirst(<(0), Iterators.map(Base.Fix1(autocor, sample), 1:length(sample))) - 1
end

function plot_autocor_time(u::Sample, k::Sample, path="tex/plots/MD_tau_energy.pdf")
    plot(; legend=:topleft, right_margin=2mm)
    xlabel!(L"$n_\textnormal{cut}$")
    ylabel!(L"$\tau$")
    for (sample, label) in
        zip((u, k, u .+ k), ("potential", "kinetic", "total") .* " energy")
        nₘₐₓ = guess_nₘₐₓ(sample)
        𝛕 = map(Base.Fix1(int_autocor_time, sample), 1:nₘₐₓ)
        plot!(1:nₘₐₓ, 𝛕; label=label)
        xlims!((1, Inf))
    end
    return savefig(path)
end
