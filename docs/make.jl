using CanopyLayers
using Documenter
using Literate
using PkgUtility




# define default docs pages
pages = Pair{Any,Any}[
    "Home" => "index.md",
    "API"  => "API.md"
]




# add example pages
gen_example = true;
gen_dir     = joinpath(@__DIR__, "src/generated");
rm(gen_dir, force=true, recursive=true);
mkpath(gen_dir);

if gen_example
    # array of example pages
    ex_pages = Pair{Any,Any}[];
    for _ex in ["bigleaf", "fluspect", "scope"]
        filename    = joinpath(@__DIR__, "src/examples/$(_ex).jl");
        script      = Literate.script(filename, gen_dir);
        code        = strip(read(script, String));
        mdpost(str) = replace(str, "@__CODE__" => code);
        Literate.markdown(filename, gen_dir, postprocess=mdpost);
        push!(ex_pages, "$(uppercasefirst(_ex))" => "generated/$(_ex).md");
    end
    # add example pages to pages
    push!(pages, "Examples" => ex_pages);
end

@info tinfo("Pages to generate:");
pretty_display!(pages);




# format the docs
mathengine = MathJax(Dict(
    :TeX => Dict(
        :equationNumbers => Dict(:autoNumber => "AMS"),
        :Macros => Dict(),
    ),
))

format = Documenter.HTML(
    prettyurls = get(ENV, "CI", nothing) == "true",
    mathengine = mathengine,
    collapselevel = 1,
    assets = ["assets/favicon.ico"]
)




# build the docs
makedocs(
    sitename = "CanopyLayers",
    format = format,
    clean = false,
    modules = [CanopyLayers],
    pages = pages,
)




# deploy the docs to Github gh-pages
deploydocs(
    repo = "github.com/Yujie-W/CanopyLayers.jl.git",
    target = "build",
    devbranch = "main",
    push_preview = true,
)
