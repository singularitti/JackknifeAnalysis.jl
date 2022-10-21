export mean, var, std, autocor, int_autocor_time, relation

mean(arg) = sum(arg) / length(arg)

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

function autocor(population::Population, n)
    μ = mean(population)
    return sum(firstindex(population):(lastindex(population) - n)) do i
        (population[i] - μ) * (population[i + n] - μ)
    end / (length(population) - n)
end

function int_autocor_time(population::Population, nₘₐₓ)
    C₀ = autocor(population, 0)
    ∑ₙ₌₁Cₙ = sum(1:nₘₐₓ) do n
        autocor(population, n)  # Cₙ
    end
    return ∑ₙ₌₁Cₙ / C₀ + 1 / 2
end

function relation(population::Population, sampler::PartitionSampler, nₘₐₓ)
    σ̂ₙ = std(population, sampler)
    τ = int_autocor_time(population, nₘₐₓ)
    s = sqrt(2τ / sampler.n) * std(Sample(population))
    return (real=σ̂ₙ, estimated=s)
end
