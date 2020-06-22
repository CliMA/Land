####
#### Defines list of tutorials given `generated_dir`
####

generate_tutorials = true

tutorials = Any[]

# Allow flag to skip generated tutorials since this is by far the slowest part of the docs build.
# generate tutorials
# following the ClimateMachine rules
if generate_tutorials
    import Literate

    include("pages_helper.jl")

    tutorials_dir = joinpath(@__DIR__, "src/tutorial_scripts")

    tutorials = [
        #=
        "Photosynthesis Testing" => [
            "Photosynthesis Model"           => "Photosynthesis/1_basics.jl",
            "Temperature Dependency"         => "Photosynthesis/2_temperature.jl",
            "Photosynthetic Rates"           => "Photosynthesis/3_photosynthesis.jl",
            "Stomatal Responses"             => "Photosynthesis/4_stomata.jl"
        ],
        =#
        "Canopy Radiative Transfer" => [
            "Reflected and Emitted Radiance" => "Radiation_Test_BRDF.jl",
            "RAMI benchmarking"              => "RAMI_benchmarking_example.jl",
            #"Soil reflectance"               => "GSV_soil_model.jl",
        ],
        "Photosynthesis" => Any[
            "Leaf Level Basics"              => "Leaf-Photosynthesis.jl",
            "Enzyme catalysis T-dependence"  => "Leaf-Photosynthesis-Rates.jl",
            "Leaf Photosynthetic rates"      => "Leaf-Photosynthesis-Synthesis.jl",
        ],
        "Canopy Energy Balance" => Any[],
    ]

    # Prepend tutorials_dir
    tutorials_jl = flatten_to_array_of_strings(get_second(tutorials))
    println("Building literate tutorials:")
    for tutorial in tutorials_jl
        println("    $(tutorial)")
    end
    transform(x) = joinpath(basename(generated_dir), replace(x, ".jl" => ".md"))
    tutorials    = transform_second(x -> transform(x), tutorials)
    tutorials_jl = map(x -> joinpath(tutorials_dir, x), tutorials_jl)

    for tutorial in tutorials_jl
        gen_dir     = joinpath(generated_dir, relpath(dirname(tutorial), tutorials_dir))
        input       = abspath(tutorial)
        script      = Literate.script(input, gen_dir)
        code        = strip(read(script, String))
        mdpost(str) = replace(str, "@__CODE__" => code)
        Literate.markdown(input, gen_dir, postprocess = mdpost)
        Literate.notebook(input, gen_dir, execute = true)
    end
end
