using JackknifeAnalysis
using LaTeXStrings
using Plots
using Plots.Measures

if !isdir("tex/plots")
    mkpath("tex/plots")
end

v1 = read("data/v1", Population)
v2 = read("data/v2", Population)
v3 = read("data/v3", Population)
v4 = read("data/v4", Population)
v5 = read("data/v5", Population)

map(mean, (v1, v2, v3, v4, v5))

let variables = (v1, v2, v3, v4, v5)
    for n in (1000, 10000)
        histograms = []
        for (i, variable) in enumerate(variables)
            X = sample(variable, ContinuousSampler(n))
            means = mean.(X)
            push!(
                histograms,
                histogram(
                    means;
                    legend=:none,
                    framestyle=:box,
                    xlims=extrema(means),
                    ylims=(0, Inf),
                    xlabel=L"$\bar{v}_%$i$",
                    ylabel=i in (1, 4) ? "frequency" : "",
                    titlefontsize=5,
                    labelfontsize=5,
                    tickfontsize=4,
                    top_margin=-1mm,
                    bottom_margin=-1mm,
                    left_margin=0mm,
                    right_margin=-1mm,
                ),
            )
            plot(histograms...)
        end
        savefig("tex/plots/sample means N=$n.pdf")
    end
end
