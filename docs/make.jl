using Documenter, GroupPresentations

DocMeta.setdocmeta!(GroupPresentations, :DocTestSetup, :(using GroupPresentations); recursive=true)

makedocs(;
    modules=[GroupPresentations],
    authors="Jean Michel <jean.michel@imj-prg.fr> and contributors",
    sitename="GroupPresentations.jl",
    format=Documenter.HTML(;
        canonical="https://jmichel7.github.io/GroupPresentations.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
    warnonly=:missing_docs,
)

deploydocs(;
    repo="github.com/jmichel7/GroupPresentations.jl",
    devbranch="main",
)
