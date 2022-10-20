function Base.read(filename::AbstractString, ::Type{Population})
    open(filename, "r") do io
        read(io, Population)
    end
end
function Base.read(io::IO, ::Type{Population})
    return Population(
        map(eachline(io)) do line
            parse(Float64, strip(line))
        end,
    )
end
