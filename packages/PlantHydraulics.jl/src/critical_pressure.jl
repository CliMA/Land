#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2022-Feb-02: migrate function to version v0.3
#     2022-Feb-02: rename the function to critical_pressure
#     2022-Feb-02: add method for LogisticVC
#     2022-Feb-02: add method for PowerVC
#     2022-Feb-02: add method for WeibullVC
#     2022-Feb-02: add method for ComplexVC
#     2022-Feb-02: add a reference kr for more customized calculations
#     2022-May-25: iterate through VCS rather than its indices for ComplexVC
#     2022-Jul-08: deflate documentations
#
#######################################################################################################################################################################################################
"""

    critical_pressure(vc::ComplexVC{FT}, kr::FT = FT(0.001)) where {FT<:AbstractFloat}
    critical_pressure(vc::LogisticVC{FT}, kr::FT = FT(0.001)) where {FT<:AbstractFloat}
    critical_pressure(vc::PowerVC{FT}, kr::FT = FT(0.001)) where {FT<:AbstractFloat}
    critical_pressure(vc::WeibullVC{FT}, kr::FT = FT(0.001)) where {FT<:AbstractFloat}

Return the critical xylem water pressure at 25 °C that triggers a given amount of loss of conductance, given
- `vc` `ComplexVC`, `LogisticVC`, `PowerVC`, or `WeibullVC` type struct
- `kr` Reference conductance, default is 0.001

"""
function critical_pressure end

critical_pressure(vc::LogisticVC{FT}, kr::FT = FT(0.001)) where {FT<:AbstractFloat} = log(kr / (vc.A + 1 - kr * vc.A)) / vc.B;

critical_pressure(vc::PowerVC{FT}, kr::FT = FT(0.001)) where {FT<:AbstractFloat} = -1 * ((1 - kr) / (kr * vc.A)) ^ (1 / vc.B);

critical_pressure(vc::WeibullVC{FT}, kr::FT = FT(0.001)) where {FT<:AbstractFloat} = -1 * (-1 * log(1 - kr)) ^ (1 / vc.C) * vc.B;

critical_pressure(vc::ComplexVC{FT}, kr::FT = FT(0.001)) where {FT<:AbstractFloat} = (
    (; VCS) = vc;

    _p_crit::FT = 0;
    for _vc in VCS
        _p_crit = min(_p_crit, critical_pressure(_vc, kr));
    end;

    return _p_crit
);
