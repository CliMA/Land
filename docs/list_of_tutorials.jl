####
#### Defines list of tutorials given `generated_dir`
####

generate_tutorials = false

tutorials = Any[]

# Allow flag to skip generated
# tutorials since this is by
# far the slowest part of the
# docs build.
if generate_tutorials

    tutorials_dir = joinpath(@__DIR__, "..", "notebooks")      # julia src files

    # generate tutorials
    import Literate

    tutorials_jl = [
        joinpath(root, f)
        for (root, dirs, files) in Base.Filesystem.walkdir(tutorials_dir)
        for f in files
    ]
    @show tutorials_jl 
    filter!(x -> endswith(x, ".jl"), tutorials_jl) # only grab .jl files
    @show tutorials_jl 
    #filter!(x -> !occursin("..jl", x), tutorials_jl)                       # currently broken, TODO: Fix me!

    println("Building literate tutorials:")
    for tutorial in tutorials_jl
        println("    $(tutorial)")
    end

    for tutorial in tutorials_jl
        gen_dir =
            joinpath(generated_dir, relpath(dirname(tutorial), tutorials_dir))
        input = abspath(tutorial)
        script = Literate.script(input, gen_dir)
        code = strip(read(script, String))
        mdpost(str) = replace(str, "@__CODE__" => code)
        Literate.markdown(input, gen_dir, postprocess = mdpost)
        # Literate.notebook(input, gen_dir, execute = true)
    end

    # TODO: Should we use AutoPages.jl?

    # These files mirror the .jl files in `CLIMA/tutorials/`:
    tutorials = Any[
        "Canopy Radiative Transfer" => Any[
            "Reflected and Emitted Radiance" => "generated/Radiation_Test_BRDF.md",
        ],
        "PhotoSynthesis" => Any[],
        "Canopy Energy Balance" => Any[],
    ]

end