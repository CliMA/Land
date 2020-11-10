using Documenter
using Literate
using WaterPhysics




# define default docs pages
pages = Any[
    "Home" => "index.md",
    "API"  => "API.md"
]




# add example pages
gen_example = true;
gen_dir     = joinpath(@__DIR__, "src/generated");
rm(gen_dir, force=true, recursive=true);
mkpath(gen_dir);

if gen_example
    filename    = joinpath(@__DIR__, "src/examples.jl");
    script      = Literate.script(filename, gen_dir);
    code        = strip(read(script, String));
    mdpost(str) = replace(str, "@__CODE__" => code);
    Literate.markdown(filename, gen_dir, postprocess=mdpost);
    push!(pages, "Examples" => "generated/examples.md");
end

@show pages;




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
)




# build the docs
makedocs(
    sitename = "WaterPhysics",
    format = format,
    clean = false,
    modules = [WaterPhysics],
    pages = pages,
)




# deploy the docs to Github gh-pages
deploydocs(
    repo = "github.com/Yujie-W/WaterPhysics.jl.git",
    target = "build",
    devbranch = "main",
    push_preview = true,
)
