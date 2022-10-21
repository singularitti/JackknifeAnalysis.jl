export Population, Sample, PartitionSampler, sample

struct Population{T} <: AbstractVector{T}
    data::Vector{T}
end
Population(A::AbstractVector) = Population(collect(A))

struct Sample{T} <: AbstractVector{T}
    data::Vector{T}
end
Sample(A::AbstractVector) = Sample(collect(A))

abstract type Sampler end
struct PartitionSampler <: Sampler
    n::Int
end

function sample(population::Population, sampler::PartitionSampler)
    return Sample.(Iterators.partition(population.data, sampler.n))
end

Base.parent(population::Population) = population.data
Base.parent(sample::Sample) = sample.data

Base.size(population::Population) = size(parent(population))
Base.size(sample::Sample) = size(parent(sample))

Base.getindex(population::Population, I...) = getindex(parent(population), I...)
Base.getindex(sample::Sample, I...) = getindex(parent(sample), I...)

Base.IndexStyle(::Type{Population{T}}) where {T} = IndexLinear()
Base.IndexStyle(::Type{Sample{T}}) where {T} = IndexLinear()
