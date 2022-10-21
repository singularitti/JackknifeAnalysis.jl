export mean, var, std, cov, cor, autocor, int_autocor_time, relation

mean(population_or_sample) = sum(population_or_sample) / length(population_or_sample)

function var(population::Population)
    μ = mean(population)
    return sum(abs2, population .- μ) / length(population)
end
function var(population::Population, sampler::PartitionSampler)
    μ = mean(population)
    samples = sample(population, sampler)
    return sum(abs2, mean(s) - μ for s in samples) / (length(samples) - 1)
end
var(sample::Sample) = var(sample, mean(sample))
var(sample::Sample, mean) = sum(abs2, sample .- mean) / (length(sample) - 1)

std(args...) = sqrt(var(args...))

function cov(population₁::Population, population₂::Population)
    if size(population₁) != size(population₂)
        throw(
            DimensionMismatch(
                "cannot compute covariance between variables of different size!"
            ),
        )
    end
    μ₁, μ₂ = mean(population₁), mean(population₂)
    return sum(zip(population₁, population₂)) do (x, y)
        (x - μ₁) * (y - μ₂)
    end / length(population₁)
end

function cor(population₁::Population, population₂::Population)
    c = cov(population₁, population₂)
    σ₁, σ₂ = std(population₁), std(population₂)
    return c / σ₁ / σ₂
end

function autocor(population::Population, n)
    μ = mean(population)
    return sum(firstindex(population):(lastindex(population) - n)) do i
        (population[i] - μ) * (population[i + n] - μ)
    end / (length(population) - n)
end

function int_autocor_time(population_or_sample, nₘₐₓ)
    C₀ = autocor(population_or_sample, 0)
    ∑ₙ₌₁Cₙ = sum(1:nₘₐₓ) do n
        autocor(population_or_sample, n)  # Cₙ
    end
    return ∑ₙ₌₁Cₙ / C₀ + 1 / 2
end

function relation(population::Population, sampler::PartitionSampler, nₘₐₓ)
    σ̂ₙ = std(population, sampler)
    τ = int_autocor_time(population, nₘₐₓ)
    s = sqrt(2τ / sampler.n) * std(Sample(population))
    return (real=σ̂ₙ, estimated=s)
end
