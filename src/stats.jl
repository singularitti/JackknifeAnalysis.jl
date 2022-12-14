export mean, var, std, cov, cor, autocor, int_autocor_time, relation

mean(data) = sum(data) / length(data)
mean(f, data) = sum(f, data) / length(data)

function var(population::Population)
    μ = mean(population)
    return sum(abs2, population .- μ) / length(population)
end
function var(population::Population, sampler::PartitionSampler)
    μ = mean(population)
    samples = sampleby(population, sampler)
    return sum(abs2, mean(s) - μ for s in samples) / (length(samples) - 1)
end
var(sample::Sample) = var(sample, mean(sample))
var(sample::Sample, mean) = sum(abs2, sample .- mean) / (length(sample) - 1)
var(sample::JackknifeSample) = var(sample, mean(sample))
function var(sample::JackknifeSample, mean)
    N = length(sample)
    return (N - 1) / N * sum(abs2, sample .- mean)
end
var(f, sample::JackknifeSample) = var(f, sample, mean(sample))
function var(f, sample::JackknifeSample, mean::Number)
    N = length(sample)
    f′ = map(f, sample)
    fₘ = f(mean)
    return (N - 1) / N * sum(abs2, f′ .- fₘ)
end
function var(f, sample::JackknifeSample, samples::JackknifeSample...)
    N = length(sample)
    if !all(length(sample′) == N for sample′ in samples)
        throw(DimensionMismatch("not all samples have the same length!"))
    end
    f′ = map(f, sample, samples...)  # f(x, y, ...) for (x, y, ...) in (sample, samples...)
    fₘ = f(map(mean, (sample, samples...))...)  # Scalar, jackknife mean = original mean
    return (N - 1) / N * sum(abs2, f′ .- fₘ)
end

std(args...) = sqrt(var(args...))

function cov(data₁, data₂)
    if size(data₁) != size(data₂)
        throw(
            DimensionMismatch("cannot compute covariance between data of different size!")
        )
    end
    μ₁, μ₂ = mean(data₁), mean(data₂)
    return sum(zip(data₁, data₂)) do (x, y)
        (x - μ₁) * (y - μ₂)
    end * _corrected(data₁)
end

function cor(data₁, data₂)
    c = cov(data₁, data₂)
    σ₁, σ₂ = std(data₁), std(data₂)
    return c / σ₁ / σ₂
end

function autocor(population::Population, n)
    μ = mean(population)
    return sum(firstindex(population):(lastindex(population) - n)) do i
        (population[i] - μ) * (population[i + n] - μ)
    end / (length(population) - n)
end
function autocor(sample::Sample, n)
    μ = mean(sample)
    return sum(firstindex(sample):(lastindex(sample) - n)) do i
        (sample[i] - μ) * (sample[i + n] - μ)
    end / (length(sample) - n - 1)
end

function int_autocor_time(data, nₘₐₓ)
    C₀ = autocor(data, 0)
    ∑ₙ₌₁Cₙ = sum(1:nₘₐₓ) do n
        autocor(data, n)  # Cₙ
    end
    return ∑ₙ₌₁Cₙ / C₀ + 1 / 2
end

function relation(population::Population, sampler::PartitionSampler, nₘₐₓ)
    σ̂ₙ = std(population, sampler)
    τ = int_autocor_time(population, nₘₐₓ)
    s = sqrt(2τ / sampler.n) * std(Sample(population))
    return (real=σ̂ₙ, estimated=s)
end
function relation(sample::Sample, nₘₐₓ)
    N = length(sample)
    σ = std(sample)
    τ = int_autocor_time(sample, nₘₐₓ)
    return sqrt(2τ / N) * σ
end

_corrected(population::Population) = inv(length(population))
_corrected(sample::Sample) = inv(length(sample) - 1)
