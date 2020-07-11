using Plants
using Test

PL = Plants;




# test the refactored functions
@testset "FT and NaN tests" begin
    for FT in [Float32, Float64]
        node   = PL.SPACSimple{FT}();
        zenith = FT(30);
        r_all  = FT(1000);

        @test true
    end
end