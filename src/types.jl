export Population, Sample, PartitionSampler, JackknifeSampler, sampleby

struct Population{T} <: AbstractVector{T}
    data::Vector{T}
end
Population(A::AbstractVector) = Population(collect(A))
function Population{S}(::UndefInitializer, n::Tuple{Int64}) where {S}
    return Population(Vector{S}(undef, n))
end

struct Sample{T} <: AbstractVector{T}
    data::Vector{T}
end
Sample(A::AbstractVector) = Sample(collect(A))
Sample{S}(::UndefInitializer, n::Tuple{Int64}) where {S} = Sample(Vector{S}(undef, n))

abstract type Sampler end
struct PartitionSampler <: Sampler
    n::Int
end
struct JackknifeSampler <: Sampler end

function sampleby(population::Population, sampler::PartitionSampler)
    return Sample.(Iterators.partition(population.data, sampler.n))
end
function sampleby(sample::Sample, ::JackknifeSampler)
    f = inv(length(sample) - 1)
    ∑ = sum(sample)
    return Sample([(∑ - value) for value in sample]) * f
end

Base.parent(population::Population) = population.data
Base.parent(sample::Sample) = sample.data

Base.size(population::Population) = size(parent(population))
Base.size(sample::Sample) = size(parent(sample))

Base.getindex(population::Population, I...) = getindex(parent(population), I...)
Base.getindex(sample::Sample, I...) = getindex(parent(sample), I...)

Base.setindex!(population::Population, v, I...) = setindex!(parent(population), v, I...)
Base.setindex!(sample::Sample, v, I...) = setindex!(parent(sample), v, I...)

Base.IndexStyle(::Type{Population{T}}) where {T} = IndexLinear()
Base.IndexStyle(::Type{Sample{T}}) where {T} = IndexLinear()

Base.similar(::Population, ::Type{S}, dims::Dims) where {S} = Population{S}(undef, dims)
Base.similar(::Sample, ::Type{S}, dims::Dims) where {S} = Sample{S}(undef, dims)
