using Documenter
using SoilHydraulics


# define default docs pages
pages = Any[
    "Home" => "index.md",
    "API"  => "API.md"
];

@show pages;


# format the docs
mathengine = MathJax(
    Dict(
        :TeX => Dict(
            :equationNumbers => Dict(:autoNumber => "AMS"),
            :Macros => Dict()
        )
    )
);

format = Documenter.HTML(
    prettyurls = get(ENV, "CI", nothing) == "true",
    mathengine = mathengine,
    collapselevel = 1,
    assets = ["assets/favicon.ico"]
);


# build the docs
makedocs(
    sitename = "SoilHydraulics",
    format = format,
    clean = false,
    modules = [SoilHydraulics],
    pages = pages
);


# deploy the docs to Github gh-pages
deploydocs(
    repo = "github.com/Yujie-W/SoilHydraulics.jl.git",
    target = "build",
    devbranch = "main",
    push_preview = true
);
