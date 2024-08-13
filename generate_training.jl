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
using CSV, DataFrames
pplt = pyimport("proplot")
pplt.close("all")
using QuasiMonteCarlo
using Surrogates
using Random

df_training = CSV.read("config/training.csv", DataFrame)

n_dim = length(df_training[:, 1])
n_samples = 1000
lb = df_training.lb
ub = df_training.ub

section("LHS Sampling")
rng = MersenneTwister(999)
_samples = sample(n_samples, lb, ub, LatinHypercubeSample(rng))
samples = zeros(n_samples, n_dim)
for i = 1:n_samples
    res = _samples[i]
    res_vec = [res[j] for j = 1:n_dim]
    samples[i, :] = res_vec
end

# samples = samples[158:159,:]
# n_samples = 2

# -17.3
# 7.5
# -19.1
# 2706.25

# -5.1
# 6.099999999999998
# 8.100000000000001
# 2293.75

# -14.299999999999999
# 18.299999999999997
# -6.099999999999998
# 1181.25

# ==============================================================================
# Settings
# ==============================================================================

# Filenames, folders
name_rotor = "djimatrice300rtk"
fname_config_rotor = "config/rotor/$name_rotor.yaml"
dir_out = "out/$name_rotor"

# Simulation conditions

rpm = samples[:, n_dim]
vRel_R = zeros(n_samples, 3) # relative velocity in rotor frame
vRel_R[:, 1] = samples[:, 1]
vRel_R[:, 2] = samples[:, 2]
vRel_R[:, 3] = samples[:, 3]
altitude = 0.0
offline_correction_rotation = false
online_correction_3D = false

# Simulation settings
lb_T = 0.0
ub_T = 250.0
tol_optim = 1e-8
maxiter_optim = 200
tol_newton = 1e-9
maxiter_newton = 1000
parallel = true


# ==============================================================================
# Setup
# ==============================================================================

config_rotor = YAML.load_file(fname_config_rotor)
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

rotor.mach_correction = false

section("Simulation")
@time res = simulate_rpm2thrust(
    rotor,
    rpm,
    vRel_R;
    lb_T=lb_T,
    ub_T=ub_T,
    parallel=parallel,
    flight_condition=flight_condition,
    tol_optim=tol_optim,
    maxiter_optim=maxiter_optim,
    trace_optim=false,
    method_optim=:brent,
    tol_newton=tol_newton,
    maxiter_newton=maxiter_newton,
);

df_training = DataFrame(
    rpm=rpm,
    u=vRel_R[:, 1],
    v=vRel_R[:, 2],
    w=vRel_R[:, 3],
    T=res.T,
    Q=res.Q,
    converged=res.converged,
)

# Convert converged to 1 and 0
df_training.converged = Int.(df_training.converged)

CSV.write(joinpath(dir_out, "trainingdata.csv"), df_training)

println(f"{sum(df_training.converged)} converged samples out of {n_samples} samples.")