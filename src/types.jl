export Population, Sample, ContinuousSampler, sample

struct Population{T} <: AbstractVector{T}
    data::Vector{T}
end
Population(A::AbstractVector) = Population(collect(A))

struct Sample{T} <: AbstractVector{T}
    data::Vector{T}
end
Sample(A::AbstractVector) = Sample(collect(A))

struct ContinuousSampler
    n::Int
end

function sample(population::Population, sampler::ContinuousSampler)
    return Sample.(Iterators.partition(population.data, sampler.n))
end

Base.parent(population::Population) = population.data
Base.parent(sample::Sample) = sample.data

Base.size(population::Population) = size(parent(population))
Base.size(sample::Sample) = size(parent(sample))

Base.getindex(population::Population, I...) = getindex(parent(population), I...)
Base.getindex(sample::Sample, I...) = getindex(parent(sample), I...)

Base.IndexStyle(::Type{<:Population}) = IndexLinear()
Base.IndexStyle(::Type{<:Sample}) = IndexLinear()
