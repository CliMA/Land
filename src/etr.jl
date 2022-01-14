"""
This function updates the electron transport rates in photosynthesis reaction center, and thus to calculate photosynthetic rate. Supported methods are

$(METHODLIST)
"""
function photosystem_electron_transport! end


"""
    photosystem_electron_transport!(ps::C3VJPModel{FT}, rc::PhotosynthesisReactionCenter{FT}, apar::FT) where {FT<:AbstractFloat}

Update the electron transport rates, given
- `ps` `C3VJPModel` type C3 photosynthesis model
- `rc` Photosynthesis reaction center
- `apar` Absorbed photosynthetically active radiation

---
# Examples
```julia
ps = C3VJPModel{Float64}();
rc = PhotosynthesisReactionCenter{Float64}();
photosystem_electron_transport!(ps, rc, 500.0);
photosystem_electron_transport!(ps, rc, 1000.0);
```
"""
photosystem_electron_transport!(ps::C3VJPModel{FT}, rc::PhotosynthesisReactionCenter{FT}, apar::FT) where {FT<:AbstractFloat} = (
    ps.j_pot = rc.F_PSII * rc.Φ_PSII_MAX * apar;
    ps.j = min(ps.j_pot, ps.j_max);

    return nothing
);


"""
    photosystem_electron_transport!(ps::C4VJPModel{FT}, rc::PhotosynthesisReactionCenter{FT}, apar::FT) where {FT<:AbstractFloat}

Update the electron transport rates, given
- `ps` `C4VJPModel` type C4 photosynthesis model
- `rc` Photosynthesis reaction center
- `apar` Absorbed photosynthetically active radiation

---
# Examples
```julia
ps = C4VJPModel{Float64}();
rc = PhotosynthesisReactionCenter{Float64}();
photosystem_electron_transport!(ps, rc, 500.0);
photosystem_electron_transport!(ps, rc, 1000.0);
```
"""
photosystem_electron_transport!(ps::C4VJPModel{FT}, rc::PhotosynthesisReactionCenter{FT}, apar::FT) where {FT<:AbstractFloat} = (
    ps.j_pot = rc.F_PSII * rc.Φ_PSII_MAX * apar;
    ps.j = ps.j_pot;

    return nothing
);
