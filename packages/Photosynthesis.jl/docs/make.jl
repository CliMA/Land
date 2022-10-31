using Documenter
using Photosynthesis


# define default docs pages
pages = Any[
    "Home" => "index.md",
    "API" => "API.md"
];


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
    sitename = "Photosynthesis",
    format = format,
    clean = false,
    modules = [Photosynthesis],
    pages = pages
);


# deploy the docs to Github gh-pages
deploydocs(
    repo = "github.com/Yujie-W/Photosynthesis.jl.git",
    target = "build",
    devbranch = "main",
    push_preview = true
);
