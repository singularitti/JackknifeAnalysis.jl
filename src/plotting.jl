using RecipesBase: @userplot, @recipe, @series

export autocorplot

@userplot AutoCorPlot
@recipe function f(plot::AutoCorPlot)
    # See http://juliaplots.org/RecipesBase.jl/stable/types/#User-Recipes-2
    sample, nₘₐₓ = plot.args
    C₀ = autocor(sample, 0)
    range = 0:nₘₐₓ
    ratios = map(range) do n
        Cₙ = autocor(sample, n)
        Cₙ / C₀
    end
    size --> (700, 400)
    markersize --> 2
    markerstrokecolor --> :auto
    markerstrokewidth --> 0
    xlims --> extrema(range)
    ylims --> (-Inf, 1)
    xguide --> raw"$n$"
    yguide --> raw"$\hat{C}(n) / \hat{C}(0)$"
    guidefontsize --> 10
    tickfontsize --> 8
    legendfontsize --> 8
    legend_foreground_color --> nothing
    legend_position --> :topright
    frame --> :box
    palette --> :tab20
    grid --> nothing
    @series begin
        seriestype --> :hline
        seriescolor --> :black
        z_order --> :back
        label --> ""
        [0]
    end
    return range, ratios
end
