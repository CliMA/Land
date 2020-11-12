using Literate




# load functions
include("pages_helper.jl")




# Allow flag to skip generated examples
generate_examples = true
if generate_examples
    examples_dir = joinpath(@__DIR__, "src/examples");
    examples     = [
        "Photosynthesis" => [
            "Temperature"    => "Photosynthesis/temperature.jl"   ,
            "Parameter Sets" => "Photosynthesis/parasets.jl"      ,
            "Photosynthesis" => "Photosynthesis/photosynthesis.jl",
        ],
        "Canopy Radiation" => [
            "Big Leaf Model" => "CanopyLayers/bigleaf.jl" ,
            "Leaf Spectrum"  => "CanopyLayers/fluspect.jl",
            "SCOPE Model"    => "CanopyLayers/scope.jl"   ,
        ],
        #=
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
        =#
    ];

    # Prepend examples_dir
    examples_jl = flatten_to_array_of_strings(get_second(examples));
    println("Building literate examples:");
    for example in examples_jl
        println("\t$(example)");
    end
    transform(x) = joinpath(basename(generated_dir), replace(x, ".jl" => ".md"));
    examples     = transform_second(x -> transform(x), examples);
    examples_jl  = map(x -> joinpath(examples_dir, x), examples_jl);

    for example in examples_jl
        gen_dir     = joinpath(generated_dir, relpath(dirname(example), examples_dir));
        input       = abspath(example);
        script      = Literate.script(input, gen_dir);
        code        = strip(read(script, String));
        mdpost(str) = replace(str, "@__CODE__" => code);
        Literate.markdown(input, gen_dir, postprocess = mdpost);
        Literate.notebook(input, gen_dir, execute = false);
    end

    # add an item into pages
    push!(pages, "Examples" => examples);
end
