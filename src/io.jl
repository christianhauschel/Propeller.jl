using .Propeller
using YAML

function get_filename_ext(path::String)
    fname = split(path, "/")[end]
    ext = split(fname, ".")[end]
    fname_noext = fname[1:end-length(ext)-1]
    return fname_noext, ext
end

"""
    export_flowunsteady(
        config_rotor::Dict, name_rotor::String, dir_out::String, rpm; 
        fc=FlightCondition(0.0), V_inf=0.0, spline_k = 5, spline_s = 1.8e-8
    )
"""
function export_flowunsteady(
    config_rotor::Dict, 
    name_rotor::String, 
    dir_out::String,
    rpm;
    fc=FlightCondition(0.0),
    V_inf=0.0,
    spline_k = 5,
    spline_s = 1.8e-8
)

    rotor = Rotor(config_rotor)
    _Re, _M, _u = estimate_Re(rotor, rpm; V_inf=V_inf, fc=fc)

    if !isdir(dir_out)
        mkpath(dir_out)
    end
    dir_out_airfoils = joinpath(dir_out, "airfoils")
    dir_out_rotors = joinpath(dir_out, "rotors")
    if !isdir(dir_out_airfoils)
        mkpath(dir_out_airfoils)
    end
    if !isdir(dir_out_rotors)
        mkpath(dir_out_rotors)
    end

    # Save config_rotor as YAML
    fname_config = joinpath(dir_out, "rotors", "$name_rotor.yaml")
    YAML.write_file(fname_config, config_rotor)


    # Rotor
    fname_blade = "$(name_rotor)_blade.csv"
    property = ["Rtip", "Rhub", "B", "blade"]
    file = [rotor.r_tip, rotor.r_hub, rotor.n_blades, fname_blade]
    description = ["(m)", "(m)", "number of blades", "blade file"]
    df_rotor = DataFrame(property=property, file=file, description=description)

    # Blade 
    property =
        ["chorddist", "pitchdist", "sweepdist", "heightdist", "airfoil_files", "spl_k", "spl_s"]
    file = [
        "$(name_rotor)_chorddist.csv",
        "$(name_rotor)_pitchdist.csv",
        "$(name_rotor)_sweepdist.csv",
        "$(name_rotor)_heightdist.csv",
        "$(name_rotor)_airfoils.csv",
        spline_k,
        spline_s,
    ]
    description = [
        "chord distribution",
        "pitch distribution",
        "LE sweep distribution",
        "LE height distribution",
        "airfoil files",
        "splines order",
        "splines smoothing",
    ]
    df_blade = DataFrame(property=property, file=file, description=description)


    α = rotor.polar_data.alpha

    _n_radial = rotor.n_radial +1

    _fnames_af = [joinpath(dir_out_airfoils, f"{replace(rotor.name, ' '=> '_')}_{i:04d}.csv") for i in 1:_n_radial]

    AF.normalize!.(rotor._airfoils)
    AF.save.(rotor._airfoils, _fnames_af)

    _fnames_af_csv_short = [get_filename_ext(fname)[1] * ".csv"  for fname in _fnames_af]

    fnames_aero = fill("", _n_radial)

    for i = 1:_n_radial
    
        # Interpolate Re 
        cl = rotor.polar_interpolator[1](α, _Re[i], rotor._s[i])
        cd = rotor.polar_interpolator[2](α, _Re[i], rotor._s[i])
        cm = rotor.polar_interpolator[3](α, _Re[i], rotor._s[i])
    
    
        fname_polar = f"s{rotor._s[i]:0.2f}_Re{_Re[i]/1e6:0.3f}.csv"
   
    
        df_aero = DataFrame("Alpha" => rad2deg.(α), "Cl" => cl, "Cd" => cd, "Cm" => cm)
        df_aero[!, "Alpha"] = (df_aero[!, "Alpha"])
    
        # if alpha[1] and alpha[end] are bigger than 90 degrees, 
        # remove all bigger than 15 degrees from the DataFrame
        if abs(df_aero[1, "Alpha"]) > 15 && df_aero[end, "Alpha"] > 15
            df_aero = df_aero[df_aero[!, "Alpha"].<15, :]
            df_aero = df_aero[df_aero[!, "Alpha"].>-15, :]
        end
    
        fnames_aero[i] = fname_polar
    
        # save df_aero as csv
        fname_aero = joinpath(dir_out_airfoils, fname_polar)
        CSV.write(fname_aero, df_aero)
    end


    println("Airfoil and polars exported to\n  $dir_out_airfoils")

    _r_nondim = rotor._r / rotor.r_tip

    # add 0 zero to the beginning of the array
    _r_nondim_af = vcat([0], _r_nondim)
    _fnames_af = vcat([_fnames_af_csv_short[1]], _fnames_af_csv_short)
    fnames_aero = vcat([fnames_aero[1]], fnames_aero)

    # round _r_nondim_af to 5 digits
    _r_nondim_af = round.(_r_nondim_af, digits=5)

    df_airfoils =
        DataFrame("r/R" => _r_nondim_af, "Contour file" => _fnames_af, "Aero file" => fnames_aero)

    _chord = FLOWMath.akima(rotor.s, rotor.chord, rotor._s)
    _theta = FLOWMath.akima(rotor.s, rotor.theta, rotor._s)
    _sweep = FLOWMath.akima(rotor.s, rotor.LEx, rotor._s)
    _height = FLOWMath.akima(rotor.s, rotor.LEz, rotor._s)

    df_chord = DataFrame("r/R" => rotor._r / rotor.r_tip, "chord/R" => _chord / rotor.r_tip)
    df_pitch = DataFrame("r/R" => rotor._r / rotor.r_tip, "twist (deg)" => rad2deg.(_theta))
    df_sweep = DataFrame(
        "r/R" => rotor._r / rotor.r_tip,
        "y/R (y-distance of LE from the middle point of hub)" => _sweep / rotor.r_tip,
    )
    df_height = DataFrame(
        "r/R" => rotor._r / rotor.r_tip,
        "z/R (height of leading edge from top face of hub)" => _height / rotor.r_tip,
    )

    # save csv
    CSV.write(joinpath(dir_out_rotors, "$name_rotor.csv"), df_rotor)
    CSV.write(joinpath(dir_out_rotors, "$(name_rotor)_blade.csv"), df_blade)
    CSV.write(joinpath(dir_out_rotors, "$(name_rotor)_chorddist.csv"), df_chord)
    CSV.write(joinpath(dir_out_rotors, "$(name_rotor)_pitchdist.csv"), df_pitch)
    CSV.write(joinpath(dir_out_rotors, "$(name_rotor)_sweepdist.csv"), df_sweep)
    CSV.write(joinpath(dir_out_rotors, "$(name_rotor)_heightdist.csv"), df_height)
    CSV.write(joinpath(dir_out_rotors, "$(name_rotor)_airfoils.csv"), df_airfoils)

    
    println("Rotor exported to\n  $dir_out_rotors")
end

