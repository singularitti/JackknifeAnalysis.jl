export Population,
    Sample, JackknifeSample, PartitionSampler, JackknifeSampler, sampleby, jackknife

abstract type Data{T<:Number} <: AbstractVector{T} end
struct Population{T} <: Data{T}
    data::Vector{T}
end
Population(A::AbstractVector) = Population(collect(A))
function Population{S}(::UndefInitializer, n::Tuple{Int64}) where {S}
    return Population(Vector{S}(undef, n))
end

struct Sample{T} <: Data{T}
    data::Vector{T}
end
Sample(A::AbstractVector) = Sample(collect(A))
Sample{S}(::UndefInitializer, n::Tuple{Int64}) where {S} = Sample(Vector{S}(undef, n))

struct JackknifeSample{T} <: Data{T}
    data::Vector{T}
end
JackknifeSample(A::AbstractVector) = JackknifeSample(collect(A))
function JackknifeSample{S}(::UndefInitializer, n::Tuple{Int64}) where {S}
    return JackknifeSample(Vector{S}(undef, n))
end

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
    return JackknifeSample([(∑ - value) for value in sample]) * f
end

jackknife(sample::Sample) = sampleby(sample, JackknifeSampler())

Base.parent(data::Data) = data.data

Base.size(data::Data) = size(parent(data))

Base.getindex(data::Data, I...) = getindex(parent(data), I...)

Base.setindex!(data::Data, v, I...) = setindex!(parent(data), v, I...)

Base.IndexStyle(::Type{<:Data}) = IndexLinear()

Base.similar(::Population, ::Type{S}, dims::Dims) where {S} = Population{S}(undef, dims)
Base.similar(::Sample, ::Type{S}, dims::Dims) where {S} = Sample{S}(undef, dims)
function Base.similar(::JackknifeSample, ::Type{S}, dims::Dims) where {S}
    return JackknifeSample{S}(undef, dims)
end

Base.BroadcastStyle(::Type{<:Data}) = Broadcast.ArrayStyle{Data}()
function Base.similar(
    bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{Data}}, ::Type{T}
) where {T}
    A = find_type(Data, bc.args)
    return similar(A, T)
end

# See https://github.com/ITensor/ITensors.jl/blob/3ac6859/src/broadcast.jl#L123-L130
function find_type(::Type{T}, args::Tuple) where {T}
    return find_type(T, find_type(T, args[1]), Base.tail(args))
end
find_type(::Type{T}, x) where {T} = x
find_type(::Type{T}, a::T, rest) where {T} = a
find_type(::Type{T}, ::Any, rest) where {T} = find_type(T, rest)
# If not found, return nothing
find_type(::Type{T}, ::Tuple{}) where {T} = nothing
