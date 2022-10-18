export readdata

function readdata(filename)
    open(filename, "r") do io
        return map(eachline(io)) do line
            parse(Float64, strip(line))
        end
    end
end
