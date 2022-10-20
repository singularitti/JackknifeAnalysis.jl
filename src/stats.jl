export mean, var, truevar, std, truestd

mean(population::Population) = sum(population) / length(population)
mean(sample::Sample) = sum(sample) / length(sample)

function var(population::Population)
    μ = mean(population)
    return sum(abs2, population .- μ) / length(population)
end
function var(sample::Sample)
    μ = mean(sample)
    return sum(abs2, sample .- μ) / (length(sample) - 1)
end

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
