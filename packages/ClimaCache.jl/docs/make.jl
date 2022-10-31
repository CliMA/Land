using ClimaCache
using Documenter


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
    sitename = "ClimaCache",
    format = format,
    clean = false,
    modules = [ClimaCache],
    pages = pages
);


# deploy the docs to Github gh-pages
deploydocs(
    repo = "github.com/Yujie-W/ClimaCache.jl.git",
    target = "build",
    devbranch = "main",
    push_preview = true
);
