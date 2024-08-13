using FlightConditions
using FLOWMath

function estimate_Re(rotor, rpm; V_inf=0.0, fc=FlightCondition(0.0))

    Ω = rpm / 60 * 2 * pi
    _u = sqrt.((Ω .* rotor._r).^2 .+ V_inf.^2)


    _Re = fc.ρ .* _u .* rotor._chord ./ fc.μ
    _M = _u ./ fc.c

    return _Re, _M, _u
end