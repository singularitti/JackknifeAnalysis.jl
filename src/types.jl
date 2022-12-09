using ConstructionBase: constructorof

export Population,
    SimpleSample, JackknifeSample, PartitionSampler, JackknifeSampler, sampleby, jackknife

# See https://github.com/JuliaLang/julia/blob/bb2d863/base/subarray.jl#L14-L24
abstract type Values{T,A<:AbstractVector{T}} <: AbstractVector{T} end
function (::Type{T})(data) where {S,A,T<:Values{S,<:A}}
    @inline
    data = map(Base.Fix1(convert, S), data)
    return constructorof(T){S,typeof(data)}(data)
end
function (::Type{T})(data) where {S,A,T<:Values{<:S,<:A}}
    @inline
    # You can't use `S` since `T` has no type parameter; therefore, `S` is not defined!
    return constructorof(T){eltype(data),typeof(data)}(data)
end
function (::Type{T})(::UndefInitializer, n) where {S,A,T<:Values{S,<:A}}
    return constructorof(T)(Vector{S}(undef, n))
end
struct Population{T,A} <: Values{T,A}
    values::A
end
abstract type Sample{T,A} <: Values{T,A} end
struct SimpleSample{T,A} <: Sample{T,A}
    values::A
end
struct JackknifeSample{T,A} <: Sample{T,A}
    values::A
end

abstract type Sampler end
struct PartitionSampler <: Sampler
    n::Int
end
struct JackknifeSampler <: Sampler end

function sampleby(population::Population, sampler::PartitionSampler)
    return SimpleSample.(Iterators.partition(population.values, sampler.n))
end
function sampleby(sample::SimpleSample, ::JackknifeSampler)
    f = inv(length(sample) - 1)
    ∑ = sum(sample)
    return JackknifeSample([(∑ - value) for value in sample]) * f
end

jackknife(sample::SimpleSample) = sampleby(sample, JackknifeSampler())

Base.parent(data::Values) = data.values

Base.size(data::Values) = size(parent(data))

Base.getindex(data::Values, I...) = getindex(parent(data), I...)

Base.setindex!(data::Values, v, I...) = setindex!(parent(data), v, I...)

Base.IndexStyle(::Type{<:Values}) = IndexLinear()

Base.similar(::Population, ::Type{S}, dims::Dims) where {S} = Population{S}(undef, dims)
Base.similar(::SimpleSample, ::Type{S}, dims::Dims) where {S} = SimpleSample{S}(undef, dims)
function Base.similar(::JackknifeSample, ::Type{S}, dims::Dims) where {S}
    return JackknifeSample{S}(undef, dims)
end

Base.BroadcastStyle(::Type{<:Values}) = Broadcast.ArrayStyle{Values}()
function Base.similar(
    bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{Values}}, ::Type{T}
) where {T}
    A = find_type(Values, bc.args)
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
