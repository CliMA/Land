"""
    get_specific_latent_heat(tem)

Specific latent heat, given
- `tem` Water temperature

Equation used is `λ = (2500.8 - 2.36*tem + 0.0016*tem^2 -0.00006*tem^3)` when tem in -25 to 40 degree C
See Polynomial curve fits to Table 2.1. R. R. Rogers; M. K. Yau (1989).
A Short Course in Cloud Physics (3rd ed.). Pergamon Press. p. 16. ISBN 0-7506-3215-1.

May need to merge with other CLIMA repository to be consistent.
"""
function get_specific_latent_heat(tem::FT) where {FT}
    temc = tem - K_0
    return slh_a0 + slh_a1*temc + slh_a2*temc^2 + slh_a3*temc^3
end
