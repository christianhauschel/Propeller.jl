"""
RPM Sweep Simulation.
"""

using FlightConditions
using Propeller
import YAML
using PyPlot, PyCall
pplt = pyimport("proplot")
pplt.close("all")

# ==============================================================================
# Settings
# ==============================================================================

# Filenames, folders
names_rotor = ["apc_845MR_model", "apc_845MR"]
fnames_config_rotor = ["config/rotor/$name_rotor.yaml" for name_rotor in names_rotor]
# dir_out = "out/$name_rotor"

# Simulation conditions
rpm = Vector(LinRange(0.0, 20000.0, 20))
n_rpm = length(rpm)
vRel_R = zeros(n_rpm, 3) # relative velocity in rotor frame

# Simulation settings
lb_T = 0.0
ub_T = 40.0
tol_optim = 1e-5
maxiter_optim = 500
parallel = true


# ==============================================================================
# Setup
# ==============================================================================

n_rotors = length(names_rotor)
rotors = Vector{Rotor}(undef, n_rotors)

for i = 1:n_rotors
    config_rotor = YAML.load_file(fnames_config_rotor[i])
    rotors[i] = Rotor(config_rotor);
    # Propeller.plot(rotors[i])
end


# ==============================================================================
# Simulation
# ==============================================================================

results = Vector{AeroResult}(undef, n_rotors)

for i = 1:n_rotors
    results[i] = simulate_rpm2thrust(
        rotors[i],
        rpm,
        vRel_R;
        lb_T = lb_T,
        ub_T = ub_T,
        tol_optim = tol_optim,
        maxiter_optim = maxiter_optim,
        parallel = parallel,
    );
end


# ==============================================================================
# Post
# ==============================================================================

fig, ax =
    pplt.subplots(figsize = (7, 3), ncols = 2, nrows = 1, sharex = true, sharey = false)
for i = 1:n_rotors
    ax[1].plot(rpm, results[i].T, "-", label=names_rotor[i])
    ax[2].plot(rpm, results[i].Q, "-", label=names_rotor[i])
end
ax[1].format(xlabel = "RPM", ylabel = "Thrust [N]", title = "Thrust vs RPM")
ax[2].format(xlabel = "RPM", ylabel = "Torque [Nm]", title = "Torque vs RPM")
ax[1].legend()
fig

# fig.savefig("out/validation/model_vs_polars.png")
