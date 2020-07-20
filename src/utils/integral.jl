###############################################################################
#
# Integral functions
#
###############################################################################
"""
    fast∫(dx::Vector{FT}, f::Vector{FT})

A fast way of integrating functions, given
- `dx` Delta x for each x
- `f` f(x) for each x
"""
function fast∫(
            dx::Vector{FT},
            f::Vector{FT}
            ) where {FT<:AbstractFloat}
    if length(dx) == length(f)
        return sum( f .* dx )
    else
        N = length(f);
        result::FT = 0.0;
        @inbounds for i=1:N
            result += f[i] * dx[i];
        end
        return result
    end
end
