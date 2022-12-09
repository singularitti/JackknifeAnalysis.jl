using JackknifeAnalysis
using Plots
using Plots.Measures

if !isdir("tex/plots")
    mkpath("tex/plots")
end

autocorplot(u, 3000; label="potential energy")
autocorplot!(k, 3000; label="kinetic energy")
autocorplot!(u .+ k, 3000; label="total energy")
savefig("tex/plots/MD_autocor_energy.pdf")

autocorplot(temperature, 300)
savefig("tex/plots/MD_autocor_temperature.pdf")

plot_autocor_v(virial, 300)
savefig("tex/plots/MD_autocor_virial.pdf")

function guess_nₘₐₓ(sample)
    return findfirst(<(0), Iterators.map(Base.Fix1(autocor, sample), 1:length(sample))) - 1
end

function plot_autocor_time_u(u::SimpleSample, path="tex/plots/MD_tau_energy.pdf")
    plot(; legend=:none, right_margin=2mm)
    xlabel!(raw"$n_\textnormal{cut}$")
    ylabel!(raw"$\tau$")
    nₘₐₓ = guess_nₘₐₓ(u)
    𝛕 = map(Base.Fix1(int_autocor_time, u), 1:nₘₐₓ)
    plot!(1:nₘₐₓ, 𝛕; label="potential energy")
    xlims!((1, Inf))
    ylims!(extrema(𝛕))
    return savefig(path)
end
function plot_autocor_time_t(
    temperature::SimpleSample, path="tex/plots/MD_tau_temperature.pdf"
)
    plot(; legend=:none, right_margin=2mm)
    xlabel!(raw"$n_\textnormal{cut}$")
    ylabel!(raw"$\tau$")
    nₘₐₓ = guess_nₘₐₓ(temperature)
    𝛕 = map(Base.Fix1(int_autocor_time, temperature), 1:nₘₐₓ)
    plot!(1:nₘₐₓ, 𝛕)
    xlims!((1, Inf))
    ylims!(extrema(𝛕))
    return savefig(path)
end
function plot_autocor_time_v(virial::SimpleSample, path="tex/plots/MD_tau_virial.pdf")
    plot(; legend=:none, right_margin=2mm)
    xlabel!(raw"$n_\textnormal{cut}$")
    ylabel!(raw"$\tau$")
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
