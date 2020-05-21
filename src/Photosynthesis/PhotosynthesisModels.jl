module PhotosynthesisModels
    using CLIMAParameters
    using DocStringExtensions
    using Parameters

    # define GLOBAL Parameters
    const GAS_R = gas_constant()
    const K_25  = 298.15 # K @ 25 Celcius

    # include the abstract types
    include("types/jmax.jl"       )
    include("types/kc.jl"         )
    include("types/ko.jl"         )
    include("types/respiration.jl")
    include("types/v_to_r.jl"     )
    include("types/vmax.jl"       )
    include("types/Γ_star.jl"     )
    include("types/para_set.jl"   )

    # include the correction equations
    include("equation/arrhenius_correction.jl"     )
    include("equation/arrhenius_peak_correction.jl")

    # include the calculation procedures
    include("procedure/get_a_v.jl"   )
    include("procedure/get_j.jl"     )
    include("procedure/get_jmax.jl"  )
    include("procedure/get_kc.jl"    )
    include("procedure/get_ko.jl"    )
    include("procedure/get_r.jl"     )
    include("procedure/get_vmax.jl"  )
    include("procedure/get_Γ_star.jl")
end