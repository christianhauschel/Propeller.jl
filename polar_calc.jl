using AirfoilPolars
import AirfoilFast as AF
using PrettySections
using PyFormattedStrings

# ==============================================================================
# Settings
# ==============================================================================

dir_data = joinpath(pwd(), "data")
dir_out = joinpath(dir_data, "polars")

af = AF.Airfoil(joinpath(dir_data, "airfoils/djimatrice300rtk/airfoil_0001.csv"))

cd_max = 1.5
n_re = 2
re = Vector(LinRange(80e3, 110e3, n_re))
alpha = Vector(-10:0.25:15)
n_re = length(re)
mach = 0.1
n_crit = 5
n_iter = 100

# ==============================================================================

polars = Vector{Polar}(undef, n_re)

init(af)

section("Solve")
for i in 1:n_re
    subsubsection(f"Re: {re[i]/1000:7.1f} k")
    polars[i] = solve(
        alpha, re[i]; 
        mach=mach, n_crit=n_crit, n_iter=n_iter,
        make_nonconverged_nan=false, interpolate_nonconverged=false,
    )
end

# polars_smooth = smooth.(polars)

fig = AirfoilPolars.plot([polars; polars_smooth])

# polars_ext = extrapolate.(polars_smooth; cd_max = cd_max )

# names_polar = generate_name.(polars_ext)

# fig = plot([polars; polars_ext]; fname="plot.png")

# fnames = [joinpath(dir_out, name*".csv") for name in names_polar]
# save.(polars_ext, fnames)

# println(fnames)

# fig