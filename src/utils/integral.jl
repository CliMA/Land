###############################################################################
#
# Integral functions
#
###############################################################################
"""
    fast∫(f::Vector{FT}, dx::Vector{FT}) where {FT<:AbstractFloat}

A fast way of integrating functions, given
- `f` f(x) for each x
- `dx` Delta x for each x
"""
function fast∫(
            f::Vector{FT},
            dx::Vector{FT}
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




"""
    fast∫!(dx::Vector{FT}, f::Vector{FT}) where {FT<:AbstractFloat}

A fast way of integrating functions, given
- `f` f(x) for each x
- `dx` Delta x for each x

Note that `f` is a local container, and its values change in this operation.
"""
function fast∫!(
            f::Array{FT,1},
            dx::Array{FT,1}
) where {FT<:AbstractFloat}
    if length(dx) == length(f)
        f .*= dx;
        return sum( f )
    else
        N = length(f);
        result::FT = 0.0;
        @inbounds for i=1:N
            result += f[i] * dx[i];
        end
        return result
    end
end
