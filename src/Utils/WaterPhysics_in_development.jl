function relative_diffusive_coefficient(T::FT) where {FT}
    return (T / K_25) ^ FT(1.8)
end

function relative_diffusive_coefficient(T::Array)
    FT = eltype(T)
    return (T ./ K_25) .^ FT(1.8)
end