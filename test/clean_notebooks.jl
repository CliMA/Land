using Test




# Use the following command to clean the notebooks if the tests fail
# jupyter nbconvert --ClearOutputPreprocessor.enabled=True --inplace *.ipynb
#
# judge if notebook is clean
@testset "Clean Notebooks" begin
    for _file in readdir(joinpath(@__DIR__, "../notebooks"))
        if _file[end-4:end] == "ipynb"
            for _line in eachline(joinpath(@__DIR__, "../notebooks", _file))
                # execution count is null?
                if occursin("   \"execution_count\":", _line)
                    @test _line == "   \"execution_count\": null,";
                end
                # metadata is empty?
                if occursin("   \"metadata\":", _line)
                    @test _line == "   \"metadata\": {},"
                end
                # output is empty?
                if occursin("   \"outputs\":", _line)
                    @test _line == "   \"outputs\": [],";
                end
            end
        end
    end
end
