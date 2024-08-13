using .Propeller

"""
    struct PolarGrid

This struct represents the aerodynamic data as gridded 3D array with the 
dimensions n_alpha × n_Re × n_span.

# Fields
- `alpha::Vector{Float64}`: Angle of attack [rad]
- `Re::Vector{Float64}`: Reynolds number [–]
- `s::Vector{String}`: Spanwise coordinate [–]
- `cl::Array{Float64, 3}`: Lift coefficient [–]
- `cd::Array{Float64, 3}`: Drag coefficient [–]
- `cm::Array{Float64, 3}`: Moment coefficient [–]
"""
struct PolarGrid
    alpha::Vector
    Re::Vector
    s::Vector
    cl::Array
    cd::Array
    cm::Array
end


function PolarGrid(filenames::Vector{Vector{String}}, span::Vector)
    n_span = length(span)
    if n_span != length(filenames)
        error("The number of polar files must be the same as the number of spanwise coordinates!")
    end

    polar = AP.load(filenames[1][1])
    alpha = deg2rad.(polar.alpha)

    n_alpha = length(alpha)
    n_Re = length(filenames[1])

    cl = Array{Float64}(undef, n_alpha, n_Re, n_span)
    cd = Array{Float64}(undef, n_alpha, n_Re, n_span)
    cm = Array{Float64}(undef, n_alpha, n_Re, n_span)
    Re = Vector{Float64}(undef, n_Re)

    __Re = Matrix{Float64}(undef, n_span, n_Re)

    # iterate files
    for i = 1:n_Re

        for j = 1:n_span
            p = AP.load(filenames[j][i])
            _Re = p.Re
            _alpha = deg2rad.(p.alpha)
            _cl = p.cl
            _cd = p.cd
            _cm = p.cm

            __Re[j, i] = _Re

            # Interpolate (as alpha is different)
            _cl = akima(_alpha, _cl, alpha)
            _cd = akima(_alpha, _cd, alpha)
            _cm = akima(_alpha, _cm, alpha)

            cl[:, i, j] = _cl
            cd[:, i, j] = _cd
            cm[:, i, j] = _cm
            Re[i] = _Re
        end

    end

    # check if re is the same along different span 
    for i = 1:n_Re
        # __Re[:,i] must be the same 
        if maximum(abs.(__Re[:, i] .- __Re[1, i])) > 1e-6
            error("Reynolds numbers must be the same along different span!")
        end
    end

    return PolarGrid(alpha, Re, span, cl, cd, cm)
end

"""
    correction3D(pg::PolarGrid, r::Vector, c::Vector, u, Ω; alpha_max_corr=30.0, alpha_linear_min=0.0, alpha_linear_max=5.0)

Apply 3D correction to the polar grid.

# Arguments
- `pg::PolarGrid`: Polar grid
- `r::Vector`: Radial coordinate [m]
- `c::Vector`: Chord length [m]
- `u::Float64`: Airfoil freestream velocity [m/s]
- `Ω::Float64`: Rotational speed [rad/s]
- `alpha_max_corr::Float64`: Maximum angle of attack for correction [deg]
- `alpha_linear_min::Float64`: Minimum angle of attack for linear correction [deg]
- `alpha_linear_max::Float64`: Maximum angle of attack for linear correction [deg]
"""
function correction3D(
    pg::PolarGrid, r::Vector, c::Vector, u::Vector, Ω;
    alpha_max_corr=30.0, alpha_linear_min=0.0, alpha_linear_max=5.0
)
    pg_3D = copy(pg)

    for i = 1:length(pg.s)
        for j = 1:length(pg.Re)
            pg_3D.cl[:, j, i], pg_3D.cd[:, j, i] = AP._correction3D(
                rad2deg.(pg.alpha), pg.cl[:, j, i], pg.cd[:, j, i], r[i], c[i], u[i], Ω;
                alpha_max_corr, alpha_linear_min, alpha_linear_max
            )
        end
    end

    return pg_3D
end

"""
    polar_interpolators(af::PolarGrid)

Return interpolators for cl, cd, and cm for a given aerodynamic data.
"""
function polar_interpolators(pg::PolarGrid)

    interpolator_cl = Interpolations.linear_interpolation(
        (pg.alpha, pg.Re, pg.s),
        pg.cl,
        extrapolation_bc=Interpolations.Flat(),
    )
    interpolator_cd = Interpolations.linear_interpolation(
        (pg.alpha, pg.Re, pg.s),
        pg.cd,
        extrapolation_bc=Interpolations.Flat(),
    )
    interpolator_cm = Interpolations.linear_interpolation(
        (pg.alpha, pg.Re, pg.s),
        pg.cm,
        extrapolation_bc=Interpolations.Flat(),
    )

    function interp_cl(alpha, Re, s)
        return interpolator_cl(alpha, Re, s)
    end

    function interp_cd(alpha, Re, s)
        return interpolator_cd(alpha, Re, s)
    end

    function interp_cm(alpha, Re, s)
        return interpolator_cm(alpha, Re, s)
    end

    return interp_cl, interp_cd, interp_cm
end

function Base.copy(af::PolarGrid)
    return PolarGrid(
        copy(af.alpha),
        copy(af.Re),
        copy(af.s),
        copy(af.cl),
        copy(af.cd),
        copy(af.cm),
    )
end


# function generate_rbg_colors(n::Int)
#     step = max(1, floor(Int, 256 / cbrt(n)))
#     colors = []

#     for r in 1:step:256
#         for g in 1:step:256
#             for b in 1:step:256
#                 push!(colors, ((r % 256 - 1) / 255, (g % 256 - 1) / 255, (b % 256 - 1) / 255))
#             end
#         end
#     end

#     return (shuffle(colors)[1:n])
# end

# function plot(p::PolarGrid; fname=nothing, dpi=300, legend=true)
#     pplt = pyimport("proplot")
#     sns = pyimport("seaborn")
#     plt = pyimport("matplotlib.pyplot")

#     fig, ax = pplt.subplots(nrows=2, ncols=2, sharex=false, sharey=false, figsize=(10, 7))

#     n_re = length(p.Re)
#     n_span = length(p.s)

#     alphas = LinRange(0.3, 1, n_re)
#     colors = generate_rbg_colors(n_span)

#     for j in 1:n_span
#         for i = 1:n_re
#             ax[1].plot(rad2deg.(p.alpha), p.cl[:, i, j], label=f"Re: {p.Re[i]/1e3:0.1f}k, s: {p.s[j]:0.2f}", lw=0.5, c=colors[j], alpha=alphas[i])
#             ax[2].plot(rad2deg.(p.alpha), p.cd[:, i, j], lw=0.5, c=colors[j], alpha=alphas[i])
#         end
#     end

#     for j in 1:n_span
#         for i = 1:n_re
#             ax[3].plot(rad2deg.(p.alpha), p.cl[:, i, j], label=f"Re = {p.Re[i]/1e3:0.1f} k, span = {p.s[j]:0.2f}", lw=0.5, c=colors[j], alpha=alphas[i])
#             ax[4].plot(rad2deg.(p.alpha), p.cd[:, i, j], lw=0.5, c=colors[j], alpha=alphas[i])
#         end
#     end

#     ax[1].set(
#         xlabel=L"$\alpha$ [deg]",
#         ylabel=L"$c_l$",
#         title="Lift Coefficient",
#     )

#     ax[2].set(
#         xlabel=L"$\alpha$ [deg]",
#         ylabel=L"$c_d$",
#         title="Drag Coefficient",
#     )

#     ax[3].set(
#         xlabel=L"$\alpha$ [deg]",
#         ylabel=L"$c_l$",
#         title="Lift Coefficient (zoom)",
#         ylim=(-1, 2),
#         xlim=(-15, 15),
#     )

#     ax[4].set(
#         xlabel=L"$\alpha$ [deg]",
#         ylabel=L"$c_d$",
#         title="Drag Coefficient (zoom)",
#         ylim = (0.0, 0.25),
#         xlim=(-15, 15),
#     )

#     if legend
#         # change legend font size 
#         ax[1].legend(ncols=2, fontsize=5)
#     end

#     if !isnothing(fname)
#         fig.savefig(fname, dpi=dpi)
#     end

#     return fig
# end
