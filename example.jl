"""
RPM Sweep Simulation.
"""

using FlightConditions
using Propeller
import YAML
using NaNStatistics
using Statistics

using PyFormattedStrings
using PrettySections
using FLOWMath

using Polynomials

# using PyPlot, PyCall
# pplt = pyimport("proplot")
# pplt.close("all")

# ==============================================================================
# Settings
# ==============================================================================

# Filenames, folders
name_rotor = "apc_11x47SF"
fname_config_rotor = "config/rotor/$name_rotor.yaml"
dir_out = "out/$name_rotor"

# Simulation conditions
n_rpm = 20
if n_rpm == 1
    rpm = Vector([2800])
else
    rpm = Vector(LinRange(1500, 4000, n_rpm))
end
vRel_R = zeros(n_rpm, 3) # relative velocity in rotor frame
vRel_R[:, 1] .= 0
vRel_R[:, 2] .= 0
vRel_R[:, 3] .= 0
altitude = 0.0
rotation_correction=:offline

# Simulation settings
lb_T = 0.0
ub_T = 100
tol_optim = 1e-8
maxiter_optim = 500
parallel = false


fit_poly = false
deg_poly = 2

plot = true


# ==============================================================================
# Setup
# ==============================================================================

config_rotor = YAML.load_file(fname_config_rotor)
config_rotor["model"]["correction"]["rotation"] = "none"
rotor = Rotor(config_rotor; rpm_offline=rpm, vRelR_offline=vRel_R);


mkpath(dir_out)

# ==============================================================================
# Simulation
# ==============================================================================

using BenchmarkTools

flight_condition = FlightCondition(altitude)
@time res = simulate_rpm2thrust(
    rotor,
    rpm,
    vRel_R;
    lb_T=lb_T,
    ub_T=ub_T,
    tol_optim=tol_optim,
);

# 20.485 ms, 11.71 MB --> none 
# 20.356 ms, 11.3 MB --> offline 
# 21.5 ms, 11.3 MB, 89.57 N, 2.9 Nm --> online

# Info 

# re_avg = zeros(n_rpm)
# α_avg = zeros(n_rpm)
# cl_avg = zeros(n_rpm)
# cd_avg = zeros(n_rpm)
# for i_rpm = 1:n_rpm
#     # global res
#     re_avg[i_rpm] = nanmean(res.re[i_rpm, :, :])
#     α_avg[i_rpm] = rad2deg.(nanmean(res.α[i_rpm, :, :]))
#     cl_avg[i_rpm] = nanmean(res.cl[i_rpm, :, :])
#     cd_avg[i_rpm] = nanmean(res.cd[i_rpm, :, :])
# end

section("Result Info")
println("RPM: ", rpm)
println("Converged: ", res.converged)
println(f"Thrust: {round.(res.T; digits=5)}")
println(f"Torque: {round.(res.Q; digits=5)}")
println(f"CT: {round.(res.CT; digits=5)}")
println(f"CQ: {round.(res.CQ; digits=5)}")
# println(f"avg. Re: {round.(re_avg; digits=1)}")
# println(f"avg. α: {round.(α_avg; digits=1)}")
# println(f"avg. cl: {round.(cl_avg; digits=3)}")
# println(f"avg. cd: {round.(cd_avg; digits=3)}")

# if fit_poly
#     poly_thrust = fit(rpm, res.T, deg_poly)
#     poly_torque = fit(rpm, res.Q, deg_poly)
# end


# ==============================================================================
# Post
# ==============================================================================

# if plot

#     fig, ax =
#         pplt.subplots(figsize=(7, 5), ncols=2, nrows=2, sharex=true, sharey=false)
#     ax[1].plot(rpm, res.T, "-")
#     ax[2].plot(rpm, res.Q, "-")
#     ax[3].plot(rpm, res.CT, "-")
#     ax[4].plot(rpm, res.CQ, "-")

#     ax[1].set(xlabel="RPM", ylabel="Thrust [N]", title="Thrust")
#     ax[2].set(ylabel="Torque [N]", title="Torque")
#     ax[3].set(ylabel=L"$C_T$ [–]", title="Thrust Coefficient")
#     ax[4].set(ylabel=L"$C_Q$ [–]", title="Torque Coefficient")

#     counter_1 = 0
#     counter_2 = 0
#     for i in 1:n_rpm
#         if res.converged[i]
#             global counter_1
#             counter_1 += 1
#             if counter_1 < 2
#                 ax[1].plot(rpm[i], res.T[i], ".", color="C0", label="Converged")
#                 ax[2].plot(rpm[i], res.Q[i], ".", color="C0", label="Converged")
#                 ax[3].plot(rpm[i], res.CT[i], ".", color="C0", label="Converged")
#                 ax[4].plot(rpm[i], res.CQ[i], ".", color="C0", label="Converged")
#             else
#                 ax[1].plot(rpm[i], res.T[i], ".", color="C0")
#                 ax[2].plot(rpm[i], res.Q[i], ".", color="C0")
#                 ax[3].plot(rpm[i], res.CT[i], ".", color="C0")
#                 ax[4].plot(rpm[i], res.CQ[i], ".", color="C0")
#             end
#         else
#             global counter_2
#             counter_2 += 1
#             if counter_2 < 2
#                 ax[1].plot(rpm[i], res.T[i], ".", color="C1", label="Not converged")
#                 ax[2].plot(rpm[i], res.Q[i], ".", color="C1", label="Not converged")
#                 ax[3].plot(rpm[i], res.CT[i], ".", color="C1", label="Not converged")
#                 ax[4].plot(rpm[i], res.CQ[i], ".", color="C1", label="Not converged")
#             else
#                 ax[1].plot(rpm[i], res.T[i], ".", color="C1")
#                 ax[2].plot(rpm[i], res.Q[i], ".", color="C1")
#                 ax[3].plot(rpm[i], res.CT[i], ".", color="C1")
#                 ax[4].plot(rpm[i], res.CQ[i], ".", color="C1")
#             end
#         end
#     end
#     ax[1].legend(ncols=1)

#     # fig.savefig("docs/img/plot.png")
#     fig.savefig(joinpath(dir_out, "plot.png"))
#     fig
# end


# Save T, Q, rpm, vRel_R to csv
# using CSV
# using DataFrames

# df = DataFrame(
#     rpm = rpm,
#     T = res.T,
#     Q = res.Q,
#     u = vRel_R[:, 1],
#     v = vRel_R[:, 2],
#     w = vRel_R[:, 3],
# )

# CSV.write(joinpath(dir_out, "data.csv"), df)