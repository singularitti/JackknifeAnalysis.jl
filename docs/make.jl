using JackknifeAnalysis
using Documenter

DocMeta.setdocmeta!(JackknifeAnalysis, :DocTestSetup, :(using JackknifeAnalysis); recursive=true)

makedocs(;
    modules=[JackknifeAnalysis],
    authors="singularitti <singularitti@outlook.com> and contributors",
    repo="https://github.com/singularitti/JackknifeAnalysis.jl/blob/{commit}{path}#{line}",
    sitename="JackknifeAnalysis.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://singularitti.github.io/JackknifeAnalysis.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/singularitti/JackknifeAnalysis.jl",
    devbranch="main",
)
