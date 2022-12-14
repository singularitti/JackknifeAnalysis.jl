using JackknifeAnalysis
using Plots

if !isdir("tex/plots")
    mkpath("tex/plots")
end

autocorplot(U, 100; label="potential energy")
autocorplot!(K, 100; label="kinetic energy")
autocorplot!(U + K, 100; label="total energy")
savefig("tex/plots/MD_autocor_energy.pdf")

autocorplot(temperature, 100; label="temperature")
savefig("tex/plots/MD_autocor_temperature.pdf")

# autocorplot(virial, 300; label="virial")
# savefig("tex/plots/MD_autocor_virial.pdf")

function guess_nₘₐₓ(sample)
    return findfirst(<(0), Iterators.map(Base.Fix1(autocor, sample), 1:length(sample))) - 1
end

autocortimeplot(U, guess_nₘₐₓ(U); label="potential energy")
savefig("tex/plots/MD_tau_energy.pdf")

autocortimeplot(temperature, guess_nₘₐₓ(temperature); label="temperature")
savefig("tex/plots/MD_tau_temperature.pdf")

# autocortimeplot(virial, guess_nₘₐₓ(virial); label="virial")
# savefig("tex/plots/MD_tau_virial.pdf")

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
