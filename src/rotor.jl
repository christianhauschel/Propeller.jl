using .Propeller
using .Propeller: PolarGrid, polar_interpolators
#using Plots

"""
    struct Rotor

This struct represents the rotor with various parameters and data for aerodynamic
analysis and simulation.

# Fields
- `name::String`: rotor name
- `r_tip::TF`: tip radius [m]
- `r_hub::TF`: hub radius [m]
- `r::Vector`: radial coordinates, mid-element [m]
- `s::Vector`: spanwise coordinates, mid-element [–]
- `_r::Vector`: radial coordinates, at node  [m]
- `_s::Vector`: spanwise coordinates, at node [–]
- `Δr::Vector`: radial segment length [m]
- `Δs::Vector`: spanwise segment length [–]
- `n_blades::Int`: number of blades [–]
- `n_radial::Int`: number of radial segments [–]
- `n_azimuthal::Int`: number of azimuthal segments [–]
- `stall_transition_rate::TF`: stall transition rate [–]
- `alpha_cutoff::TF`: stall cut-off AoA [deg]
- `alpha_cl0::TF`: zero-lift AoA [deg]
- `theta::Vector`: twist angle, mid-element [deg]
- `_theta::Vector`: twist angle, at node [rad]
- `chord::Vector`: chord length, mid-element [m]
- `_chord::Vector`: chord length, at node [m]
- `LEx::Vector`: FLOWUnsteady sweep (x-distance of LE from the middle point of hub), mid-element [m]
- `_LEx::Vector`: FLOWUnsteady sweep (x-distance of LE from the middle point of hub), at node [m]
- `LEz::Vector`: FLOWUnsteady dihedral (height of leading edge from top face of hub), mid-element [m]
- `_LEz::Vector`: FLOWUnsteady dihedral (height of leading edge from top face of hub), at node [m]
- `cl_a_slope::Vector`: lift slope [1/rad]
- `r_nondim::Vector`: nondimensional radial coordinates [–]
- `psi::Vector`: azimuthal coordinates [rad]
- `polar_data::PolarGrid`: polar data
- `polar_interpolator::Tuple{Function, Function}`: polar interpolators
- `use_polar_model::Bool`: use polar model
- `_airfoils::Vector{AF.Airfoil}`: airfoils at nodes
- `area_crosssection::Vector`: airfoil area, mid-element [m²]
"""
mutable struct Rotor
    name::String
    r_tip
    r_hub
    r::Vector
    s::Vector
    _r::Vector
    _s::Vector
    Δr::Vector
    Δs::Vector
    n_blades
    n_radial
    n_azimuthal
    stall_transition_rate # stall transition rate
    alpha_cutoff # stall cut-off AoA
    alpha_cl0 # zero-lift AoA
    theta::Vector
    _theta::Vector
    chord::Vector
    _chord::Vector
    LEx::Vector
    _LEx::Vector
    LEz::Vector 
    _LEz::Vector
    cl_a_slope::Vector
    r_nondim::Vector
    psi::Vector

    polar_data::Union{PolarGrid,Nothing}
    polar_interpolator::Union{Tuple{Function,Function,Function},Nothing}

    online_correction3D::Bool
    mach_correction::Bool 
    hubtip_correction::Symbol 
    inflow_model::Symbol

    use_polar_model::Bool
    _airfoils::Vector{AF.Airfoil}
    area_crosssection::Vector
end

function Base.copy(rotor::Rotor)
    return Rotor(
        copy(rotor.name),
        copy(rotor.r_tip),
        copy(rotor.r_hub),
        copy(rotor.r),
        copy(rotor.s),
        copy(rotor._r),
        copy(rotor._s),
        copy(rotor.Δr),
        copy(rotor.Δs),
        copy(rotor.n_blades),
        copy(rotor.n_radial),
        copy(rotor.n_azimuthal),
        copy(rotor.stall_transition_rate),
        copy(rotor.alpha_cutoff),
        copy(rotor.alpha_cl0),
        copy(rotor.theta),
        copy(rotor._theta),
        copy(rotor.chord),
        copy(rotor._chord),
        copy(rotor.LEx),
        copy(rotor._LEx),
        copy(rotor.LEz),
        copy(rotor._LEz),
        copy(rotor.cl_a_slope),
        copy(rotor.r_nondim),
        copy(rotor.psi),
        copy(rotor.polar_data),

        copy(rotor.polar_interpolator),
        copy(rotor.use_polar_model),

        copy(rotor.online_correction3D),
        copy(rotor.mach_correction),
        copy(rotor.hubtip_correction),
        copy(rotor.inflow_model),

        copy(rotor._airfoils),
        copy(rotor.area_crosssection),
    )
end


function Base.show(io::IO, rotor::Rotor)
    println(io, f"──── Rotor ────")
    println(io, f" name:       {rotor.name}")
    println(io, f" r_tip:       {rotor.r_tip:5.2f}")
    println(io, f" r_hub:       {rotor.r_hub:5.2f}")
    println(io, f" n_blades:    {rotor.n_blades:5d}")
    println(io, f" n_radial:    {rotor.n_radial:5d}")
    println(io, f" n_azimuthal: {rotor.n_azimuthal:5d}")
end

function update_chord!(rotor::Rotor, s::Vector, chord::Vector)
    rotor.chord = akima(s, chord, rotor.s)
    rotor._chord = akima(s, chord, rotor._s)
end

function update_twist!(rotor::Rotor, s::Vector, theta::Vector)
    rotor.theta = akima(s, theta, rotor.s)
    rotor._theta = akima(s, theta, rotor._s)
end

"""
    Rotor(config::Dict; dummy=false, rotation_correction=:none, rpm_offline=nothing)

This function creates a `Rotor` object based on the configuration dictionary.

# Arguments
- `config::Dict`: configuration dictionary
- `dummy::Bool`: dummy rotor
- `correction3D::Symbol`: :none, :online, :offline
"""
function Rotor(config::Dict; dummy=false, rpm_offline=nothing, vRelR_offline=nothing)
    rotation_correction = Symbol(config["model"]["correction"]["rotation"])

    hubtip_correction = Symbol(config["model"]["correction"]["hubtip"])
    mach_correction = config["model"]["correction"]["mach"]
    inflow_model = Symbol(config["model"]["inflow"])

    n_radial = config["n_radial"]
    n_azimuthal = config["n_azimuthal"]
    n_blades = config["n_blades"]
    r_tip = config["r_tip"]
    r_hub = config["r_hub"]
    _cl_a_slope = config["polars"]["model"]["cl_a_slope"]
    use_polar_model = config["polars"]["model"]["enable"]
    stall_transition_rate = float(config["polars"]["model"]["stall_transition_rate"])
    alpha_cutoff = float(config["polars"]["model"]["alpha_cutoff"])
    alpha_cl0 = float(config["polars"]["model"]["alpha_cl0"])

    airfoils_raw = AF.Airfoil.(config["airfoils"]["fnames"])

    # Radial Sampling
    if "radial_sampling" in keys(config)
        sampling = Symbol(config["radial_sampling"]["type"])
        if "coeff_conical" in keys(config["radial_sampling"])
            coeff_conical_sampling = config["radial_sampling"]["coeff_conical"]
        else
            coeff_conical_sampling = 1
        end
        if "ratio_geometric" in keys(config["radial_sampling"])
            shrink_ratio_geometric = config["radial_sampling"]["ratio_geometric"]
        else
            shrink_ratio_geometric = 1 / 1.2
        end
    else
        sampling = :uniform
        shrink_ratio_geometric = 1 / 1.2
        coeff_conical_sampling = 1
    end
    r, s, Δr, Δs, _r, _s = radial_segments(
        n_radial, r_hub, r_tip;
        sampling=sampling,
        coeff_conical_sampling=coeff_conical_sampling,
        shrink_ratio_geometric=shrink_ratio_geometric,
    )

    r_nondim = r ./ r_tip

    chord = akima(config["chord"]["span"], config["chord"]["values"], s)

    # FLOWUnsteady sweep/dihedral
    if "LEx" in keys(config)
        LEx = akima(config["LEx"]["span"], config["LEx"]["values"], s)
    else 
        LEx = zeros(length(s))
    end 
    if "LEz" in keys(config)
        LEz = akima(config["LEz"]["span"], config["LEz"]["values"], s)
    else 
        LEz = zeros(length(s))
    end 

    if dummy == false
        _polar_data = PolarGrid(config["polars"]["fnames"], config["polars"]["span"])
        _polar_interpolator = polar_interpolators(_polar_data)

        if rotation_correction == :offline || rotation_correction == :online
            cl = zeros(length(_polar_data.alpha), length(_polar_data.Re), n_radial)
            cd = zeros(length(_polar_data.alpha), length(_polar_data.Re), n_radial)
            cm = zeros(length(_polar_data.alpha), length(_polar_data.Re), n_radial)

            for i = 1:length(_polar_data.alpha)
                for j = 1:length(_polar_data.Re)
                    for k = 1:n_radial
                        cl[i, j, k] = _polar_interpolator[1](_polar_data.alpha[i], _polar_data.Re[j], s[k])
                        cd[i, j, k] = _polar_interpolator[2](_polar_data.alpha[i], _polar_data.Re[j], s[k])
                        cm[i, j, k] = _polar_interpolator[3](_polar_data.alpha[i], _polar_data.Re[j], s[k])
                    end
                end
            end
            polar_data = PolarGrid(_polar_data.alpha, _polar_data.Re, s, cl, cd, cm)
            polar_interpolator = polar_interpolators(polar_data)

            if rotation_correction == :offline
                if rpm_offline === nothing || vRelR_offline === nothing
                    @warn "No RPM or relative velocity was specified for offline rotation correction. Skipping correction!"
                else
                    r_polar = akima(s,r, _polar_data.s)
                    chord_polar = akima(s, chord, _polar_data.s)

                    Ω = mean(rpm_offline) / 60 * 2π
                    # Approximate airfoil freestream velocity 
                    Vwind = mean(-vRelR_offline[3])
                    if sign(Vwind) == -1
                        Vwind = 0
                    end
                    u = (Vwind .^ 2 .+ (Ω .* r_polar) .^ 2) .^ 0.5

                    _polar_data = correction3D(_polar_data, r_polar, chord_polar, u, Ω)
                end
            end
        else 
            polar_data = _polar_data
            polar_interpolator = _polar_interpolator
        end
    else
        polar_data = nothing
        polar_interpolator = nothing
    end

    psi = Vector(LinRange(0, 2π, n_azimuthal))


    _chord = akima(s, chord, _s)
    theta = akima(config["twist"]["span"], config["twist"]["values"], s)
    _theta = akima(s, theta, _s)
    _LEx = akima(s, LEx, _s)
    _LEz = akima(s, LEz, _s)

    cl_a_slope = akima(
        config["polars"]["model"]["cl_a_slope"]["span"],
        config["polars"]["model"]["cl_a_slope"]["values"],
        s,
    )

    # Interpolate airfoils & Calculate area
    if length(airfoils_raw) != length(config["airfoils"]["span"])
        error("Number of airfoils does not match the number of spanwise coordinates!")
    end
    _airfoils = AF.interpolate_airfoils(airfoils_raw, config["airfoils"]["span"], _s)
    airfoils = AF.interpolate_airfoils(airfoils_raw, config["airfoils"]["span"], s)

    area_crosssection = AF.area.(airfoils) .* chord .^ 2

    return Rotor(
        config["name"],
        r_tip,
        r_hub,
        r,
        s,
        _r,
        _s,
        Δr,
        Δs,
        n_blades,
        n_radial,
        n_azimuthal,
        stall_transition_rate, # stall transition rate
        alpha_cutoff, # stall cut-off AoA
        alpha_cl0, # zero-lift AoA
        theta,
        _theta,
        chord,
        _chord,
        LEx,
        _LEx,
        LEz,
        _LEz,
        cl_a_slope,
        r_nondim,
        psi,

        polar_data,
        polar_interpolator,

        rotation_correction==:online,
        mach_correction,
        hubtip_correction,
        inflow_model,

        use_polar_model,
        _airfoils,
        area_crosssection,
    )
end



function r_hub(rotor::Rotor)
    return rotor.r_nondim[1] * rotor.r_tip
end



"""
This function calculates the spanwise coordinate based on the radial coordinate.
"""
function span2radius(span, r_shroud, hubtip_ratio)
    return span * r_shroud * (1 - hubtip_ratio) + r_shroud * hubtip_ratio
end

"""
This function calculates the spanwise coordinate based on the radial coordinate.
"""
function radius2span(r, r_hub, r_tip)
    return (r - r_hub) / (r_tip - r_hub)
end


"""
Calculate the midpoints of a vector.

The resulting vector has one element less than the input vector.
"""
function midpoints(x)
    return 0.5 .* (x[1:end-1] + x[2:end])
end


"""
    radial_segments(n_radial::Int, r_hub, r_tip)

This function calculates the radial and spanwise coordinates of the rotor. The
radial coordinates are calculated at the mid-element of each segment. The spanwise
coordinates are calculated at the mid-element of each segment. The radial and
spanwise coordinates are also calculated at the nodes of each segment.

# Arguments
- `n_radial::Int`: number of radial segments [–]
- `r_hub`: hub radius [m]
- `r_tip`: tip radius [m]
- `sampling::Symbol`: radial sampling method
- `coeff_conical_sampling::Float64`: coefficient for conical sampling
- `shrink_ratio_geometric::Float64`: shrinkage ratio for geometric sampling


# Returns
- `r::Vector`: radial coordinates, mid-element [m]
- `s::Vector`: spanwise coordinates, mid-element [–]
- `_r::Vector`: radial coordinates, at node  [m]
- `_s::Vector`: spanwise coordinates, at node [–]
- `Δr::Vector`: radial segment length [m]
- `Δs::Vector`: spanwise segment length [–]
"""
function radial_segments(
    n_radial, r_hub, r_tip;
    sampling=:uniform,
    coeff_conical_sampling=1,
    shrink_ratio_geometric=1 / 1.2
)::Tuple
    # Radius
    if sampling == :uniform
        _r = Vector(LinRange(r_hub, r_tip, n_radial + 1))
    elseif sampling == :cosine
        _r = cosine(r_hub, r_tip, n_radial + 1)
    elseif sampling == :linear
        _r = conical(r_hub, r_tip, n_radial + 1; coeff=coeff_conical_sampling)
    elseif sampling == :geometric
        _r = geometric(r_hub, r_tip, n_radial + 1; q=shrink_ratio_geometric)
    elseif sampling == :geometric2
        _r = geometric2(r_hub, r_tip, n_radial + 1; r=shrink_ratio_geometric)
    else
        error("Invalid radial distribution!")
    end

    r = midpoints(_r)

    # Span
    _s = radius2span.(_r, r_hub, r_tip)
    s = radius2span.(r, r_hub, r_tip)

    # Calculate segments' length
    Δr = diff(_r)
    Δs = diff(_s)

    return r, s, Δr, Δs, _r, _s
end

"""
    subdivide_span(_s::Vector, n_between::Int)::Vector

Subdivide the spanwise coordinates by adding `n_between` more points in-between.

# Arguments
- `_s::Vector`: span at nodes [–]
- `n_between_sections::Int`: number of points to add in-between
"""
function subdivide_span(_s::Vector, n_between::Int)::Vector
    _n_radial_old = length(_s)
    n_radial_old = _n_radial_old - 1

    _n_radial = n_between * n_radial_old + 1
    n_radial = _n_radial - 1

    # Calculate new spanwise coordinates
    _s_new = zeros(_n_radial)
    k = 2
    for i = 2:n_between:_n_radial
        _s_new[i] = (_s[k] + _s[k-1]) / 2
        k += 1
    end
    k = 1
    for i = 1:n_between:_n_radial
        _s_new[i] = _s[k]
        k += 1
    end
    return _s_new
end

"""
    linear_subdivide(nodes::Vector{T}, n::Int) where T

Subdivide the nodes by adding `n` more points in-between, linearly.
"""
function linear_subdivide(nodes::Vector{T}, n::Int) where {T}
    result = Vector{T}(undef, (length(nodes) - 1) * n + 1)

    for i in 1:length(nodes)-1
        delta = (nodes[i+1] - nodes[i]) / n
        result[(i-1)*n+1] = nodes[i]

        for j in 2:n
            result[(i-1)*n+j] = nodes[i] + delta * (j - 1)
        end
    end

    result[end] = nodes[end]

    return result
end

function plot(rotor::Rotor; fname=nothing)
    p = Plots.plot(layout=(2, 2), size=(700, 600), legend=:none)
    plot!(p[1], rotor._r, rotor._chord, marker = :circle, markerstrokewidth = 0, xlabel="radius [m]", ylabel="chord [m]", label=rotor.name)
    # plot!(p[2], rotor._r/rotor.r_tip, rotor._chord./rotor.r_tip, marker = :circle, markerstrokewidth = 0, xlabel="radius/R", ylabel="chord/R")
    plot!(p[2], rotor._s, rad2deg.(rotor._theta), marker = :circle, markerstrokewidth = 0, xlabel="span", ylabel="twist [deg]")
    plot!(p[3], rotor._s, rotor._LEx, marker = :circle, markerstrokewidth = 0, xlabel="span", ylabel="LEx [m]")
    plot!(p[4], rotor._s, rotor._LEz, marker = :circle, markerstrokewidth = 0, xlabel="span", ylabel="LEz [m]")
    plot!(p[1], legend=:topright)

    if fname !== nothing
        savefig(p, fname)
    end
    return p
end

# function plot(rotors::Vector{Rotor}; fname=nothing, dpi=300)
    
#     p = Plots.plot(layout=(2, 2), size=(700, 600), legend=:none, dpi=dpi)
#     for rotor in rotors
#         plot!(p[1], rotor._r, rotor._chord, marker = :circle, markerstrokewidth = 0, xlabel="radius [m]", ylabel="chord [m]", label=rotor.name)
#         # plot!(p[2], rotor._r/rotor.r_tip, rotor._chord./rotor.r_tip, marker = :circle, markerstrokewidth = 0, xlabel="radius/R", ylabel="chord/R")
#         plot!(p[2], rotor._s, rad2deg.(rotor._theta), marker = :circle, markerstrokewidth = 0, xlabel="span", ylabel="twist [deg]")
#         plot!(p[3], rotor._s, rotor._LEx, marker = :circle, markerstrokewidth = 0, xlabel="span", ylabel="LEx [m]")
#         plot!(p[4], rotor._s, rotor._LEz, marker = :circle, markerstrokewidth = 0, xlabel="span", ylabel="LEz [m]")
#     end

#     # legend to p[1]
#     plot!(p[1], legend=:topright)

#     if fname !== nothing
#         savefig(p, fname)
#     end
#     return p
# end

# ==============================================================================
# Corrections
# ==============================================================================


# """
#     correction3D(polar_data::PolarGrid, r_nondim::Vector, chord::Vector, tsr::Float64)

# This function applies 3D correction to the polar data. The correction is based on
# the tip speed ratio (TSR) ratio. 

# The TSR is calculated based on the mean RPM and axial velocity (averaged over time).
# """
# function correction3D!(polar_data, rpm::Vector, vRel_R::Matrix)
#     Ω = mean(rpm) / 60 * 2π

#     # Approximate airfoil freestream velocity 
#     Vwind = mean(-vRel_R[3])
#     if sign(Vwind) == -1
#         Vwind = 0
#     end
#     u = (Vwind.^2 .+ (Ω .* rotor.r).^2).^0.5

#     polars_new = correction3D(polar_data, rotor.r, rotor.chord, u, Ω)

#     plot([polar_data[end], polars_new[end]])

#     # rotor.polar_data = polars_new
# end


