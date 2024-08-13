"""
RPM Sweep Simulation.
"""

using FlightConditions
using Propeller
import YAML
using NaNStatistics
using Statistics
using PyPlot, PyCall
using PyFormattedStrings
using PrettySections
using Polynomials
using XLSX, DataFrames

pplt = pyimport("proplot")
pplt.close("all")


# ==============================================================================
# Settings
# ==============================================================================

# Filenames, folders
# name_rotor = "djimatrice300rtk"
name_rotor = "apc_11x47SF"
fname_config_rotor = "config/rotor/$name_rotor.yaml"
dir_out = "out/$name_rotor"

# Simulation conditions
n_rpm = 20
rpm = Vector(LinRange(1000, 10000, n_rpm))
vRel_R = zeros(n_rpm, 3) # relative velocity in rotor frame
vRel_R[:, 1] .= 0.0
vRel_R[:, 2] .= 0.0
vRel_R[:, 3] .= 0.0
altitude = 0.0
offline_correction_rotation = false
online_correction_3D = true

# Simulation settings
lb_T = 0.0
ub_T = 150.0
tol_optim = 1e-8
maxiter_optim = 1000
parallel = true


fit_poly = true
deg_poly = 2

plot = true

# Validation data
xlsx_data = XLSX.readtable("config/rotor/apc_11x47SF/performance.xlsx", "ref_data")
df = DataFrame(xlsx_data)
df.RPM = convert(Array{Float64}, df.RPM)
df.T = convert(Array{Float64}, df.T)
df.Q = convert(Array{Float64}, df.Q)
df.CT = convert(Array{Float64}, df.CT)
df.V = convert(Array{Float64}, df.V)

# convert df.V to float 
df = df[df.V .== 0.0, :]
df = df[df.RPM .<= maximum(rpm), :]


# ==============================================================================
# Setup
# ==============================================================================

config_rotor = YAML.load_file(fname_config_rotor)
config_rotor["model"]["correction"]["rotation"] = "none"
config_rotor["model"]["correction"]["mach"] = true
rotor = Rotor(config_rotor);


if offline_correction_rotation
    println("Apply rotation correction to rotor polars.")
    correction3D!(rotor, rpm, vRel_R)
end

mkpath(dir_out)

# ==============================================================================
# Simulation
# ==============================================================================

flight_condition = FlightCondition(altitude)

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



# ==============================================================================
# Post
# ==============================================================================

if plot

    fig, ax =
        pplt.subplots(figsize=(7, 5), ncols=2, nrows=2, sharex=true, sharey=false)
    ax[1].plot(rpm, res.T, "-", color="C0")
    ax[2].plot(rpm, res.Q, "-", color="C0")
    ax[3].plot(rpm, res.CT, "-", color="C0")
    ax[4].plot(rpm, res.CQ, "-", color="C0")

    

    ax[1].set(xlabel="RPM", ylabel="Thrust [N]", title="Thrust")
    ax[2].set(ylabel="Torque [N]", title="Torque")
    ax[3].set(ylabel=L"$C_T$ [–]", title="Thrust Coefficient")
    ax[4].set(ylabel=L"$C_Q$ [–]", title="Torque Coefficient")

    counter_1 = 0
    counter_2 = 0
    for i in 1:n_rpm
        if res.converged[i]
            global counter_1
            counter_1 += 1
            if counter_1 < 2
                ax[1].plot(rpm[i], res.T[i], ".", color="C0", label="Converged")
                ax[2].plot(rpm[i], res.Q[i], ".", color="C0", label="Converged")
                ax[3].plot(rpm[i], res.CT[i], ".", color="C0", label="Converged")
                ax[4].plot(rpm[i], res.CQ[i], ".", color="C0", label="Converged")
            else
                ax[1].plot(rpm[i], res.T[i], ".", color="C0")
                ax[2].plot(rpm[i], res.Q[i], ".", color="C0")
                ax[3].plot(rpm[i], res.CT[i], ".", color="C0")
                ax[4].plot(rpm[i], res.CQ[i], ".", color="C0")
            end
        else
            global counter_2
            counter_2 += 1
            if counter_2 < 2
                ax[1].plot(rpm[i], res.T[i], ".", color="C1", label="Not converged")
                ax[2].plot(rpm[i], res.Q[i], ".", color="C1", label="Not converged")
                ax[3].plot(rpm[i], res.CT[i], ".", color="C1", label="Not converged")
                ax[4].plot(rpm[i], res.CQ[i], ".", color="C1", label="Not converged")
            else
                ax[1].plot(rpm[i], res.T[i], ".", color="C1")
                ax[2].plot(rpm[i], res.Q[i], ".", color="C1")
                ax[3].plot(rpm[i], res.CT[i], ".", color="C1")
                ax[4].plot(rpm[i], res.CQ[i], ".", color="C1")
            end
        end
    end

    ax[1].plot(df.RPM, df.T, "--", color="C1", label="Ref")
    ax[2].plot(df.RPM, df.Q, "--", color="C1", label="Ref")
    ax[3].plot(df.RPM, df.CT, "--", color="C1", label="Ref")

    ax[1].legend(ncols=1)

    fig.savefig(joinpath("docs/img", "apc11x47.png"));
    fig
end


