using AirfoilPolars

# ==============================================================================
# Settings
# ==============================================================================

dir_data = "/mnt/a/Code/10_drones/Propeller/data"
dir_polars = joinpath(dir_data, "polars")
# name_af = "NACA4408"
# fname = joinpath(dir_polars, f"{name_af}_Re0.020_M0.00_N9.0.csv")
fname = joinpath(dir_polars, "NACA_4412_Re0.010_M0.10_N7.0_360_V.plr")
cd_max = 1.5

# ==============================================================================

p = load(fname)
# p.name_airfoil = name_af
p_smooth = smooth(p)
p_ext = extrapolate(p_smooth; cd_max=cd_max)

fig = plot([p_ext, p])

name = generate_name(p_ext)
save(p_ext, joinpath(dir_polars, name*".csv"))

fig