"""
RPM Sweep Simulation.
"""

using FlightConditions
using Propeller
import YAML
using NaNStatistics
using Statistics
using XLSX, DataFrames
using PyFormattedStrings
using PrettySections
using FLOWMath

using Polynomials

using PyPlot, PyCall
pplt = pyimport("proplot")
pplt.close("all")
plt = pyimport("matplotlib.pyplot")
np = pyimport("numpy")

# ==============================================================================
# Settings
# ==============================================================================

# Filenames, folders
# name_rotor = "djimatrice300rtk"
name_rotor = "dji9443"
# name_rotor = "rotor_anopp"
# name_rotor = "apc_11x47SF"
# name_rotor = "rotor_simple"
fname_config_rotor = "config/rotor/$name_rotor.yaml"
dir_out = "out/$name_rotor"

rpm_base = 5400
Ω_base = rpm_base * 2 * pi / 60
period_base = 1 / (rpm_base / 60)
n_revs = 1

fs = 500
Δt = 1 / fs
tend = period_base * n_revs
Δangle_deg = rad2deg(Δt * Ω_base)

t = collect(0:Δt:tend)
n_t = length(t)

# Simulation conditions

if n_t == 1
    rpm = Vector([rpm_base])
else
    rpm = Vector(LinRange(rpm_base, rpm_base, n_t))
end
vRel_R = zeros(n_t, 3) # relative velocity in rotor frame
vRel_R[:, 1] .= 0.0
vRel_R[:, 2] .= 0
vRel_R[:, 3] .= 0.0
# vRel_R[:,3] .= 0.0
altitude = 0.0


# Simulation settings
lb_T = 0.0
ub_T = 10
tol_optim = 1e-8
maxiter_optim = 500
parallel = true

mach_correction = true

fit_poly = false
deg_poly = 2

plot = true


# ==============================================================================
# Setup
# ==============================================================================

config_rotor = YAML.load_file(fname_config_rotor)
config_rotor["n_azimuthal"] = 120
config_rotor["n_radial"] = 30
rotor = Rotor(config_rotor);


mkpath(dir_out)

# ==============================================================================
# Simulation
# ==============================================================================

flight_condition = FlightCondition(89726.67666252334, 291.64)
# flight_condition = FlightCondition(altitude)

using BenchmarkTools



@time res = simulate_rpm2thrust(
    rotor,
    rpm,
    vRel_R;
    lb_T=lb_T,
    ub_T=ub_T,
    tol_optim=tol_optim,
    maxiter_optim=maxiter_optim,
    parallel=parallel,
    flight_condition=flight_condition,
);

# grid, _, _ = blade2grid(rotor; compact=true)
# @time Np, Tp = map_loads(rotor, rpm, t, res, grid; angle0=0.0);


# fig, ax = pplt.subplots(figsize=(7, 5), ncols=1, nrows=2, sharex=true, sharey=false)
# # ax[1].plot(t, Np[:, 1, 1], "-")
# # ax[1].plot(t, Np[:, 2, 1], "-")
# # ax[1].plot(t, Np[:, 5, 1], "-")
# ax[1].plot(t, Np[:, 11, 1], "-")
# # ax[2].plot(t, Tp[:, 1, 1], "-")
# # ax[2].plot(t, Tp[:, 2, 1], "-")
# # ax[2].plot(t, Tp[:, 5, 1], "-")
# ax[2].plot(t, Tp[:, 11, 1], "-")
# ax[1].set(
#     xlabel="t [s]",
#     ylabel="Np [N/m]",
# )
# ax[2].set(
#     xlabel="t [s]",
#     ylabel="Tp [N/m]",
# )
# fig

# df_val_Np = DataFrame(XLSX.readtable("validation/loading/loading.xlsx", "Np"))
# df_val_Tp = DataFrame(XLSX.readtable("validation/loading/loading.xlsx", "Tp"))

fig, ax = pplt.subplots(
    figsize=(6, 4),
    nrows=2,
    sharex=true, sharey=false
)
ax[1].plot(rotor.r/rotor.r_tip, res.Np[1,:,1], label="Hybrid BEM")
# ax[1].plot(df_val_Np[!, "r/R"], df_val_Np[!, "Np"], label="Validation (FLOW unsteady)")
ax[1].set(
    xlabel="r/R",
    ylabel="Np [N/m]",
)

ax[2].plot(rotor.r/rotor.r_tip, res.Tp[1,:,1], label="Hybrid BEM")
# ax[2].plot(df_val_Tp[!, "r/R"], df_val_Tp[!, "Tp"], label="Validation (FLOW unsteady)")
ax[2].set(
    xlabel="r/R",
    ylabel="Tp [N/m]",
)

ax[1].legend()
# fig.savefig("docs/img/loads.png", dpi=300)
fig