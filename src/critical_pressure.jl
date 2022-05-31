#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2022-Feb-02: migrate function to version v0.3
#     2022-Feb-02: rename the function to critical_pressure
#
#######################################################################################################################################################################################################
"""
This function returns the critical xylem water pressure (at 25 °C) that triggers 99.triggers a given amount of loss of hydraulic conductance. The supported methods are

$(METHODLIST)

"""
function critical_pressure end


#######################################################################################################################################################################################################
#
# Changes to the method
# General
#     2022-Feb-02: add method for ComplexVC
#     2022-Feb-02: add a reference kr for more customized calculations
# Bug fixes:
#     2022-May-25: iterate through VCS than than its indices
#
#######################################################################################################################################################################################################
"""

    critical_pressure(vc::ComplexVC{FT}, kr::FT = FT(0.001)) where {FT<:AbstractFloat}

Return the critical xylem water pressure at 25 °C that triggers a given amount of loss of conductance, given
- `vc` `ClimaCache.ComplexVC` type struct
- `kr` Reference conductance, default is 0.001
"""
critical_pressure(vc::ComplexVC{FT}, kr::FT = FT(0.001)) where {FT<:AbstractFloat} = (
    @unpack VCS = vc;

    _p_crit::FT = 0;
    for _vc in VCS
        _p_crit = min(_p_crit, critical_pressure(_vc, kr));
    end;

    return _p_crit
);


#######################################################################################################################################################################################################
#
# Changes to the method
# General
#     2022-Feb-02: add method for LogisticVC
#     2022-Feb-02: add a reference kr for more customized calculations
#
#######################################################################################################################################################################################################
"""

    critical_pressure(vc::LogisticVC{FT}, kr::FT = FT(0.001)) where {FT<:AbstractFloat}

Return the critical xylem water pressure at 25 °C that triggers a given amount of loss of conductance, given
- `vc` `ClimaCache.LogisticVC` type struct
- `kr` Reference conductance, default is 0.001
"""
critical_pressure(vc::LogisticVC{FT}, kr::FT = FT(0.001)) where {FT<:AbstractFloat} = (
    @unpack A, B = vc;

    return log(kr / (A + 1 - kr * A)) / B
);


#######################################################################################################################################################################################################
#
# Changes to the method
# General
#     2022-Feb-02: add method for PowerVC
#     2022-Feb-02: add a reference kr for more customized calculations
#
#######################################################################################################################################################################################################
"""

    critical_pressure(vc::PowerVC{FT}, kr::FT = FT(0.001)) where {FT<:AbstractFloat}

Return the critical xylem water pressure at 25 °C that triggers a given amount of loss of conductance, given
- `vc` `ClimaCache.PowerVC` type struct
- `kr` Reference conductance, default is 0.001
"""
critical_pressure(vc::PowerVC{FT}, kr::FT = FT(0.001)) where {FT<:AbstractFloat} = (
    @unpack A, B = vc;

    return -1 * ((1 - kr) / (kr * A)) ^ (1 / B)
);


#######################################################################################################################################################################################################
#
# Changes to the method
# General
#     2022-Feb-02: add method for WeibullVC
#     2022-Feb-02: add a reference kr for more customized calculations
#
#######################################################################################################################################################################################################
"""

critical_pressure(vc::WeibullVC{FT}, kr::FT = FT(0.001)) where {FT<:AbstractFloat}

Return the critical xylem water pressure at 25 °C that triggers a given amount of loss of conductance, given
- `vc` `ClimaCache.WeibullVC` type struct
- `kr` Reference conductance, default is 0.001
"""
critical_pressure(vc::WeibullVC{FT}, kr::FT = FT(0.001)) where {FT<:AbstractFloat} = (
    @unpack B, C = vc;

    return -1 * (-1 * log(1 - kr)) ^ (1 / C) * B
);
