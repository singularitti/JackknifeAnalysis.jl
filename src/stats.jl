export mean, var, truevar, std, truestd, autocor

mean(population::Population) = sum(population) / length(population)
mean(sample::Sample) = sum(sample) / length(sample)

function var(population::Population)
    μ = mean(population)
    return sum(abs2, population .- μ) / length(population)
end
var(sample::Sample) = var(sample, mean(sample))
var(sample::Sample, mean) = sum(abs2, sample .- mean) / (length(sample) - 1)

function truevar(population::Population, sampler::ContinuousSampler)
    μ = mean(population)
    samples = sample(population, sampler)
    return sum(abs2, mean(s) - μ for s in samples) / (length(samples) - 1)
end

std(population::Population) = sqrt(var(population))
std(sample::Sample) = sqrt(var(sample))

function truestd(population::Population, sampler::ContinuousSampler)
    return sqrt(truevar(population, sampler))
end

function autocor(population::Population, n)
    μ = mean(population)
    return sum(firstindex(population):(lastindex(population) - n)) do i
        (population[i] - μ) * (population[i + n] - μ)
    end / (length(population) - n)
end
