using RecipesBase: @userplot, @recipe, @series

export autocorplot, autocortimeplot

@userplot AutoCorPlot
@recipe function f(plot::AutoCorPlot)
    # See http://juliaplots.org/RecipesBase.jl/stable/types/#User-Recipes-2
    sample, nₘₐₓ = plot.args
    C₀ = autocor(sample, 0)
    terms = 0:nₘₐₓ
    ratios = map(terms) do n
        Cₙ = autocor(sample, n)
        Cₙ / C₀
    end
    size --> (700, 400)
    markersize --> 2
    markerstrokecolor --> :auto
    markerstrokewidth --> 0
    xlims --> extrema(terms)
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
        seriestype --> :scatter
        terms, ratios
    end
    @series begin
        seriestype --> :path
        z_order --> :back
        label := ""
        terms, ratios
    end
end

@userplot AutocorTimePlot
@recipe function f(plot::AutocorTimePlot)
    sample, nₘₐₓ = plot.args
    terms = 1:nₘₐₓ
    𝛕 = map(Base.Fix1(int_autocor_time, sample), terms)
    size --> (700, 400)
    markersize --> 2
    markerstrokecolor --> :auto
    markerstrokewidth --> 0
    xlims --> extrema(terms)
    ylims --> extrema(𝛕)
    xguide --> raw"$n_\textnormal{cut}$"
    yguide --> raw"integrated autocorrelation time ($\tau$)"
    guidefontsize --> 10
    tickfontsize --> 8
    legendfontsize --> 8
    legend_foreground_color --> nothing
    legend_position --> :bottomright
    frame --> :box
    palette --> :tab20
    grid --> nothing
    @series begin
        seriestype --> :scatter
        terms, 𝛕
    end
    @series begin
        seriestype --> :path
        z_order --> :back
        label := ""
        terms, 𝛕
    end
end
