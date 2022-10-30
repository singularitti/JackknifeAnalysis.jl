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

function plot_autocor_u(
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
function plot_autocor_t(
    temperature::Sample, nₘₐₓ=3000, path="tex/plots/MD_autocor_temperature.pdf"
)
    plot(; legend=:none, right_margin=3mm)
    xlabel!(L"$n$")
    ylabel!(L"$\hat{C}(n) / \hat{C}(0)$")
    C₀ = autocor(temperature, 0)
    N = 0:nₘₐₓ
    r = map(N) do n
        Cₙ = autocor(temperature, n)
        Cₙ / C₀
    end
    plot!(N, r)
    xlims!(extrema(N))
    return savefig(path)
end
function plot_autocor_v(virial::Sample, nₘₐₓ=3000, path="tex/plots/MD_autocor_virial.pdf")
    plot(; legend=:none, right_margin=3mm)
    xlabel!(L"$n$")
    ylabel!(L"$\hat{C}(n) / \hat{C}(0)$")
    C₀ = autocor(virial, 0)
    N = 0:nₘₐₓ
    r = map(N) do n
        Cₙ = autocor(virial, n)
        Cₙ / C₀
    end
    plot!(N, r)
    xlims!(extrema(N))
    return savefig(path)
end

function guess_nₘₐₓ(sample)
    return findfirst(<(0), Iterators.map(Base.Fix1(autocor, sample), 1:length(sample))) - 1
end

function plot_autocor_time_u(u::Sample, path="tex/plots/MD_tau_energy.pdf")
    plot(; legend=:none, right_margin=2mm)
    xlabel!(L"$n_\textnormal{cut}$")
    ylabel!(L"$\tau$")
    nₘₐₓ = guess_nₘₐₓ(u)
    𝛕 = map(Base.Fix1(int_autocor_time, u), 1:nₘₐₓ)
    plot!(1:nₘₐₓ, 𝛕; label="potential energy")
    xlims!((1, Inf))
    ylims!(extrema(𝛕))
    return savefig(path)
end
function plot_autocor_time_t(temperature::Sample, path="tex/plots/MD_tau_temperature.pdf")
    plot(; legend=:none, right_margin=2mm)
    xlabel!(L"$n_\textnormal{cut}$")
    ylabel!(L"$\tau$")
    nₘₐₓ = guess_nₘₐₓ(temperature)
    𝛕 = map(Base.Fix1(int_autocor_time, temperature), 1:nₘₐₓ)
    plot!(1:nₘₐₓ, 𝛕)
    xlims!((1, Inf))
    ylims!(extrema(𝛕))
    return savefig(path)
end
function plot_autocor_time_v(virial::Sample, path="tex/plots/MD_tau_virial.pdf")
    plot(; legend=:none, right_margin=2mm)
    xlabel!(L"$n_\textnormal{cut}$")
    ylabel!(L"$\tau$")
    nₘₐₓ = guess_nₘₐₓ(virial)
    𝛕 = map(Base.Fix1(int_autocor_time, virial), 1:nₘₐₓ)
    plot!(1:nₘₐₓ, 𝛕; label="virial")
    xlims!((1, Inf))
    ylims!(extrema(𝛕))
    return savefig(path)
end

function covariance_matrix(variables)
    m = [cov(variable₁, variable₂) for variable₁ in variables, variable₂ in variables]
    return map(x -> round(x; digits=5), m)
end

function correlation_matrix(variables)
    ρ = [cor(variable₁, variable₂) for variable₁ in variables, variable₂ in variables]
    return map(x -> round(x; digits=5), ρ)
end

getbins(sample, b=19) = Sample(mean.(sampleby(Population(sample), PartitionSampler(b))))

pressure(T, v, Nₚ=1000) = 1 + inv(3 * Nₚ * T) * v

function get_f_std(f, variable₁, variable₂, b=19)
    sample₁, sample₂ = map(
        variable -> jackknife(getbins(variable, b)), (variable₁, variable₂)
    )
    return std(f, sample₁, sample₂)
end