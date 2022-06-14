#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jun-13: add function to interpolate the spectrum
#
#######################################################################################################################################################################################################
"""
This function interpolate the spectrum to give values at the target wavelength bin(s). The supported methods include

$(METHODLIST)

"""
function read_spectrum end


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jun-13: add method to interpolate the spectrum
#
#######################################################################################################################################################################################################
"""

    read_spectrum(x::Vector{FT}, y::Vector{FT}, target::FT) where {FT<:AbstractFloat}

Return the spectrum value at target wavelength bin, given
- `x` X-axis of the spectrum
- `y` Y-axis of the spectrum
- `target` Target x value
"""
read_spectrum(x::Vector{FT}, y::Vector{FT}, target::FT) where {FT<:AbstractFloat} = (
    @assert length(x) == length(y) "Dimensions of provided spectrum x and y must match!";
    @assert x[1] <= target <= x[end] "Target wavelength must be within the range provided spectum!";

    # iterate through the spectrum and find the index
    _ind = 0;
    for _i in 1:length(x)-1
        if x[_i] <= target <= x[_i+1]
            _ind = _i;
            break;
        end;
    end;

    return ((x[_ind+1] - target) * y[_ind] + (target - x[_ind]) * y[_ind+1]) / (x[_ind+1] - x[_ind])
);


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jun-13: add method to interpolate the spectrum via multiple steps
#
#######################################################################################################################################################################################################
"""

    read_spectrum(x::Vector{FT}, y::Vector{FT}, x₁::FT, x₂::FT; steps::Int = 2) where {FT<:AbstractFloat}

Return the spectrum value at target wavelength bin, given
- `x` X-axis of the spectrum
- `y` Y-axis of the spectrum
- `x₁` Lower x boundary
- `x₂` Upper x boundary
- `steps` The incremental Δx is `(x₂ - x₁) / steps`
"""
read_spectrum(x::Vector{FT}, y::Vector{FT}, x₁::FT, x₂::FT; steps::Int = 2) where {FT<:AbstractFloat} = (
    _xs = collect(FT,range(x₁, x₂; length=steps+1));
    _ys = read_spectrum.([x], [y], _xs);

    return mean(_ys)
);
