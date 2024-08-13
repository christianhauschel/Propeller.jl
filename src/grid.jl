using .Propeller

function blade2grid(
    rotor::Rotor;
    chordwise_position=0.0,
    n_between_sections=1,
    ratio_expansion=30.0,
    n_upper=5,
    n_lower=5,
    x_direction=1,
    angle0=0.0,
    compact=false,
)

    if n_between_sections > 1
        error("not implemented yet!")
    end

    _n_radial = rotor.n_radial + 1


    # Discretization between each airfoil
    sect = [
        [(1.0, n_between_sections, 1.0, false)] for i = 1:rotor.n_radial
    ]
    expansion_ratio = ratio_expansion * ones(_n_radial)

    # Orientation of chord of each airfoil (yaw, pitch, roll)
    orien = zeros((rotor.n_radial + 1, 3))
    orien[:, 1] = rad2deg.(rotor._theta)

    n_pts = n_upper + n_lower + 2

    _airfoils = AF.scale.(rotor._airfoils, rotor._chord)

    _origin = zeros((_n_radial, 3))
    if typeof(chordwise_position) == String
        if chordwise_position == "centroid"
            centers = AF.centroid.(_airfoils)
            C = hcat(centers...)'

            _origin[:, 1] = C[:, 1]
            _origin[:, 3] = C[:, 2]
        end
    else
        _origin[:, 1] = chordwise_position .* rotor._chord
    end
    _origin[:, 2] = rotor._r

    grid = nothing
    x_over_c_upper = nothing
    x_over_c_lower = nothing




    if !compact
        # ======================================================================
        # Lofted Blade
        # ======================================================================

        grid = zeros(3, _n_radial, n_pts, 1)
        x_over_c_upper = zeros(_n_radial, n_upper + 1)
        x_over_c_lower = zeros(_n_radial, n_lower + 1)

        for i in 1:_n_radial
            x = (_airfoils[i].x .- _origin[i, 1]) * x_direction
            y = _airfoils[i].y .- _origin[i, 3]

            # Separate upper and lower sides to make the contour injective in x
            upper, lower = gt.splitcontour(x, y)

            # Parameterize both sides independently
            fun_upper =
                gt.parameterize(upper[1], upper[2], zeros(Float64, size(upper[1])); inj_var=1)
            fun_lower =
                gt.parameterize(lower[1], lower[2], zeros(Float64, size(lower[1])); inj_var=1)

            # New discretization for both surfaces
            upper_points = gt.discretize(fun_upper, 0, 1, n_upper, expansion_ratio[i]; central=true)
            lower_points = gt.discretize(fun_lower, 0, 1, n_lower, expansion_ratio[i]; central=true)



            # Put both surfaces back together from TE over the top and from LE over the bottom.
            reverse!(upper_points)                           # Trailing edge over the top
            new_x_upper = [point[1] for point in upper_points]
            new_y_upper = [point[2] for point in upper_points]    # Leading edge over the bottom
            new_x_lower = [point[1] for point in lower_points]
            new_y_lower = [point[2] for point in lower_points]

            chord_upper = maximum(new_x_upper) - minimum(new_x_upper)
            chord_lower = maximum(new_x_lower) - minimum(new_x_lower)
            x_over_c_upper[i, :] = new_x_upper ./ chord_upper
            x_over_c_lower[i, :] = new_x_lower ./ chord_lower


            new_x = vcat(new_x_upper, new_x_lower)
            new_y = vcat(new_y_upper, new_y_lower)

            airfoil_pts = zeros(n_pts, 3)
            for k = 1:n_pts
                airfoil_pts[k, :] = [new_x[k], rotor._r[i], new_y[k]]

                # Positions the airfoil along the blade in the right orientation
                rot_matrix = gt.rotation_matrix(0, x_direction * rad2deg(rotor._theta[i]), 0)
                airfoil_pts[k, :] = gt.transform(airfoil_pts[k, :], rot_matrix, [0, 0, 0])

                # Rotate around z axis 
                rot_matrix = gt.rotation_matrix(rad2deg(angle0), 0, 0)
                airfoil_pts[k, :] = gt.transform(airfoil_pts[k, :], rot_matrix, [0, 0, 0])
            end
            grid[:, i, :, 1] = airfoil_pts'
        end


    else
        # ======================================================================
        # Compact Blade
        # ======================================================================

        grid = zeros(3, _n_radial, 1, 1)

        lifting_line = zeros(_n_radial, 3)
        lifting_line[:, 1] = (0.25 * rotor._chord .- _origin[:, 1]) * x_direction
        lifting_line[:, 2] = _origin[:, 2]
        lifting_line[:, 3] = 0.0 .- _origin[:, 3]

        rot_matrix_z = gt.rotation_matrix(rad2deg(angle0), 0, 0)
        for i = 1:_n_radial
            # Rotate around y axis
            rot_matrix_y = gt.rotation_matrix(0, x_direction * rad2deg(rotor._theta[i]), 0)
            lifting_line[i, :] = gt.transform(lifting_line[i, :], rot_matrix_y, [0, 0, 0])

            # Rotate around z axis 
            lifting_line[i, :] = gt.transform(lifting_line[i, :], rot_matrix_z, [0, 0, 0])
        end

        grid[:, :, 1, 1] = lifting_line'
        x_over_c_upper = zeros(_n_radial, 1)
        x_over_c_lower = zeros(_n_radial, 1)
    end

    return grid, x_over_c_upper, x_over_c_lower
end


function map_loads(rotor::Rotor, rpm::Vector, t::Vector, res::AeroResult, grid::Array{Float64,4}; angle0=0.0)

    # Dimensions of the grid
    dims = size(grid)
    # n_t = dims[1]
    n_t = length(t)
    _n_radial = dims[2]
    n_chord = dims[3]

    compact = n_chord == 1

    if !compact
        error("not implemented yet!")
    end


    # Rotate grid 
    Ω = rpm / 60.0 * 2.0 * pi
    angle = cumul_integrate(t, Ω) .+ angle0

    grid_temporal = zeros(n_t, 3, _n_radial, n_chord)
    for i = 1:n_t
        R = gt.rotation_matrix(rad2deg(angle[i]), 0, 0)
        for j = 1:_n_radial
            for k = 1:n_chord
                grid_temporal[i, :, j, k] = R * grid[:, j, k, 1]
            end
        end
    end

    # Calculate angles of grid points
    angles = zeros(n_t, _n_radial, n_chord)
    for i = 1:n_t
        for j = 1:_n_radial
            for k = 1:n_chord
                _angles = atan(-grid_temporal[i, 1, j, k], grid_temporal[i, 2, j, k]) # TODO: check signs here

                if _angles < 0
                    _angles += 2π
                end

                angles[i, j, k] = _angles
            end
        end
    end

    # fig, ax = pplt.subplots()
    # ax[1].plot(t, rad2deg.(angles[:, 11, 1]), "-")
    # fig.savefig("angles.png")
    # fig


    

    # Interpolate loads to grid
    Np = zeros(n_t, _n_radial, n_chord)
    Tp = zeros(n_t, _n_radial, n_chord)
    n_azimuthal = size(res.Np)[3]
    θ = LinRange(0, 2π, n_azimuthal)
    for i = 1:n_t

        # Extrapolate loading to radial grid points
        _Np = zeros(_n_radial, n_azimuthal)
        _Tp = zeros(_n_radial, n_azimuthal)
        for l = 1:n_azimuthal
            interpNp = linear_interpolation(
                rotor.s, res.Np[i, :, l],
                extrapolation_bc = Flat(),
            )
            interpTp = linear_interpolation(
                rotor.s, res.Tp[i, :, l],
                extrapolation_bc = Flat(),
            )
            _Np[:, l] = interpNp(rotor._s)
            _Tp[:, l] = interpTp(rotor._s)
        end

        for j = 1:_n_radial
            for k = 1:n_chord
                # Interpolate loading to azimuthal grid points
                Np[i, j, k] = linear(θ, _Np[j, :], angles[i, j, k])
                Tp[i, j, k] = linear(θ, _Tp[j, :], angles[i, j, k])
            end
        end
    end

    return Np, Tp

end

