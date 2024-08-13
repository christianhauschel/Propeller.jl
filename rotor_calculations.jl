using Propeller
import YAML
using PrettySections
using DataFrames
using FLOWMath
using FlightConditions
import AirfoilFast as AF

name_rotor = "dji_9443" # 40k ... 300k
fname_config_rotor = "config/rotor/$name_rotor.yaml"

# ==============================================================================
# Estimate Reynolds number
# ==============================================================================

fc = FlightCondition(0.0)
rpm = 3000
dummy = true # set to true if polars/airfoils not defined

config_rotor = YAML.load_file(fname_config_rotor)
rotor = Rotor(config_rotor; dummy=dummy)
V_inf = 10.0

Re, M, u = estimate_Re(rotor, rpm; V_inf=V_inf, fc=fc)


# Create a DataFrames
df = DataFrame(s = rotor._s, r = rotor._r, Re = Re, u = u, chord = rotor._chord, ρ = fc.ρ, μ = fc.μ, M=M)

section("Results")
println(df)


# radius at 75% span
chord75 = akima(rotor.s, rotor.chord, 0.75)

ratio_chord75_R = chord75 / rotor.r_tip


_t_max = AF.thickness_max.(rotor._airfoils) .* rotor._chord

s_ref = config_rotor["twist"]["span"]

t_max_ref = akima(rotor._s, _t_max, s_ref)