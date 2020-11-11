###############################################################################
#
# Recursive FT test
#
###############################################################################
"""
    FT_test(para::Any, FT)

Test the the floating point type of para is FT, given
- `para` Any type of parameters
- `FT` given FT

If data type is not supported, use pass!
"""
function FT_test(para::Array, FT)
    passed = true;

    # fail if para is float but not FT
    if eltype(para) <: AbstractFloat
        if eltype(para) != FT
            passed = false;
        end
    else
        if !all(FT_test.(para, FT))
            passed = false;
        end
    end

    return passed
end




function FT_test(para::DataType)
    return true
end




function FT_test(para::Function)
    return true
end




function FT_test(para::Module)
    return true
end




function FT_test(para::Number, FT)
    passed = true;

    # fail if para is float but not FT
    if typeof(para) <: AbstractFloat
        if typeof(para) != FT
            passed = false;
        end
    end

    return passed
end




function FT_test(para::Symbol)
    return true
end




function FT_test(para::Any, FT)
    passed = true;

    # try to detech struct
    try
        arr = [];
        for fn in fieldnames( typeof(para) )
            push!(arr, FT_test( getfield(para, fn), FT ));
        end
        if !all(arr)
            passed = false;
        end
    catch e
        nothing
    end

    return passed
end








###############################################################################
#
# Recursive NaN test
#
###############################################################################
"""
    NaN_test(para::Any)

Test the the floating point type of para is not NaN, given
- `para` Any type of parameters

If data type is not supported, use pass!
"""
function NaN_test(para::Array)
    passed = true;

    # fail if para is number
    if eltype(para) <: Number
        if !all(.!isnan.(para))
            passed = false;
        end
    else
        if !all(NaN_test.(para))
            passed = false;
        end
    end

    return passed
end




function NaN_test(para::DataType)
    return true
end




function NaN_test(para::Function)
    return true
end




function NaN_test(para::Module)
    return true
end




function NaN_test(para::Number)
    passed = true;

    # fail if para is NaN
    if isnan(para)
        passed = false;
    end

    return passed
end




function NaN_test(para::Symbol)
    return true
end




function NaN_test(para::Any)
    passed = true;

    # try to detech struct
    try
        arr = [];
        for fn in fieldnames( typeof(para) )
            push!(arr, NaN_test( getfield(para, fn) ));
        end
        if !all(arr)
            passed = false;
        end
    catch e
        nothing
    end

    return passed
end
