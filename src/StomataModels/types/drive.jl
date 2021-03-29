###############################################################################
#
# Driver for leaf gas exchange
#
###############################################################################
"""
    abstract type AbstractDrive

Hierarchy of AbstractDrive
- [`GlcDrive`](@ref)
- [`GswDrive`](@ref)
"""
abstract type AbstractDrive end




"""
    struct GlcDrive

Gas exchange update is driven by changes in total leaf diffusive conductance to
    CO₂
"""
struct GlcDrive <: AbstractDrive end




"""
    struct GswDrive

Gas exchange update is driven by changes in stomtal conductance to H₂O
"""
struct GswDrive <: AbstractDrive end
