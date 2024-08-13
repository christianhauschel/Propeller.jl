"""
Comparison with data from Fig. 9 in [1].

# References
[1] B. Davoudi, “A Hybrid Blade Element Momentum Model for Flight Simulation of 
    Rotary Wing Unmanned Aerial Vehicles,” AIAA Paper, 2019.
"""

using Revise
using PrettySections
using FlightConditions
using Propeller
using PyCall, PyPlot
using YAML
using PyFormattedStrings    
pplt = pyimport("proplot")
pplt.close("all")

fname_config = "config/rotor/apc_845MR.yaml"


# flight_condition = FlightCondition(101325, 293.15)
flight_condition = FlightCondition(0.0)


n_rpm = 30
α_R = [0, 30, 60, 89.9]
lb_T = 0.0
ub_T = 10.0
# α_R = [89.9]
V_inf = [15]
sheetname = f"APC 8x4.5MR_omega_{V_inf[1]}"
mach_correction = false

rpm = Vector(LinRange(0, 10000, n_rpm))
n_V = length(V_inf)
n_α_R = length(α_R)

section("Setup")

config_rotor = YAML.load_file(fname_config)
subsubsection("Rotor")
@time "TOTAL" rotor = Rotor(config_rotor);

# Propeller.plot(rotor)

T = zeros((n_α_R, n_V, n_rpm))
V_rel_B = zeros((n_α_R, n_V, n_rpm, 3))
μ = zeros((n_α_R, n_V, n_rpm))

converged = zeros(Bool, (n_α_R, n_V, n_rpm))

section("Simulation")

rotor.mach_correction = true
# rotor.hubtip_correction = :none

res = nothing

@time for i_α = 1:n_α_R
    for i_V = 1:n_V
        V_rel_B[i_α, i_V, :, 1] .= -V_inf[i_V] * cosd(α_R[i_α])
        # V_rel_B[i_α, i_V, :, 2] .= -V_inf[i_V] * cosd(α_R[i_α])
        V_rel_B[i_α, i_V, :, 3] .= -V_inf[i_V] * sind(α_R[i_α])
        global res
        res = simulate_rpm2thrust(
            rotor,
            rpm,
            V_rel_B[i_α, i_V, :,:];
            lb_T = lb_T,
            ub_T = ub_T,
            flight_condition = flight_condition,
        )

        T[i_α, i_V, :] = res.T  
        converged[i_α, i_V, :] = res.converged
        
        for i_rpm = 1:n_rpm
            if converged[i_α, i_V, i_rpm] == false
                T[i_α, i_V, i_rpm] = NaN
            end
        end
    end
end

# ==============================================================================
# Plots
# ==============================================================================

section("Post")

using DataFrames
using XLSX
df = DataFrame(XLSX.readtable("validation/davoudi/validation.xlsx", sheetname))


fig, ax = pplt.subplots(figsize = (6, 4))
for i_α = 1:n_α_R
    for i_V = 1:n_V
        my_filter(alpha, V) = alpha == α_R[i_α] && V == abs.(V_inf[i_V])
        ax[1].plot(rpm, T[i_α, i_V, :], "-", label = f"Sim ({α_R[i_α]:0.1f}°, {V_inf[i_V]:0.1f} m/s)", color = "C$i_α")
        _df = filter([:alpha_R, :V_inf] => my_filter, df)
        ax[1].plot(_df.rpm, _df.T, "o", color = "C$i_α", label=f"Exp ({α_R[i_α]:0.1f}°, {V_inf[i_V]:0.1f} m/s)")
    end
end
ax[1].legend(ncols = 1)
ax[1].set(title = "Validation", xlabel = "RPM [1/min]", ylabel = "T [N]")
fig.savefig(f"out/validation_davoudi/validation_V{V_inf[1]}.png")
fig
