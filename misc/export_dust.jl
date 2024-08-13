"""Convert Propeller.Rotor to Dust"""


using FlightConditions
using Propeller
import YAML
import FLOWMath as FM
using NaNStatistics
using Statistics
using Trapz
using PyPlot, PyCall
using PyFormattedStrings
using PrettySections
import AirfoilFast as AF
pplt = pyimport("proplot")
pplt.close("all")

# ==============================================================================
# Settings
# ==============================================================================

# Filenames, folders
name_rotor = "apc_845MR_model"
fname_config_rotor = "config/rotor/$name_rotor.yaml"
dir_out = "out/$name_rotor"
dir_af = joinpath(dir_out, "airfoils")
dir_af_short = "airfoils"
fname_blade = joinpath(dir_out, "blade.in")

# Dust settings
n_span = 1
nelem_chord = 40
type_chord = "cosineLE"
type_span = "uniform"
reference_chord_fraction = 0.0 # 0.5 = half of the chord, centroid: mean centroid of the airfoils
meshfiletype = "parametric"
elType = "p"


# ==============================================================================
# Load Rotor
# ==============================================================================

config_rotor = YAML.load_file(fname_config_rotor)
rotor = Rotor(config_rotor);


# ==============================================================================
# Write File
# ==============================================================================

# make dir_af if not existent
if !isdir(dir_af)
    mkpath(dir_af)
end


_n_radial = length(rotor._r)

_chord = rotor._chord 
_twist = rad2deg.(rotor._theta)


_sweep = zeros(_n_radial)
_dihed = zeros(_n_radial)

function dirname(fname)
    return joinpath(split(fname, "/")[1:end-1]...)
end

_fnames_af = [joinpath(dir_af, f"{i:02d}.dat") for i in 1:_n_radial]
_fnames_af_short = [joinpath(dir_af_short, f"{i:02d}.dat") for i in 1:_n_radial]

_airfoils = rotor._airfoils
AF.normalize!.(_airfoils)
AF.save_dust.(_airfoils, _fnames_af)

_r = rotor._r

if reference_chord_fraction == "centroid"
    _centroids = AF.centroid.(rotor._airfoils)
    _centroids = hcat(_centroids...)'
    chord_position = Trapz.trapz(rotor._r, _centroids[:, 1]) / (rotor.r_tip - rotor.r_hub)
else
    chord_position = reference_chord_fraction
end

open(fname_blade, "w") do f
    write(f, f"meshFileType = {meshfiletype}\n")
    write(f, f"elType = {elType}\n")
    write(f, "\n")
    write(f, "starting_point = (/0.0, $(rotor.r_hub), 0.0/) !Beginning of the extruded geometry\n")
    write(f, f"reference_chord_fraction = {chord_position}\n")
    write(f, f"nelem_chord = {nelem_chord}\n")
    write(f, "\n")
    write(f, f"type_chord = {type_chord}\n")
    write(f, "\n\n")
    write(f, "!-----------------------------------------\n")
    write(f, "! Sections\n")
    write(f, "!-----------------------------------------\n\n")

    for i in 1:_n_radial
        write(f, "! Section $(i)\n")
        write(f, "chord = $(_chord[i])\n")
        write(f, "twist = $(_twist[i])\n")
        write(f, "airfoil = $(_fnames_af_short[i])\n")
        write(f, "\n")

        if i < _n_radial
            write(f, "! Region $(i)–$(i+1)\n")
            write(f, "span = $(_r[i+1] - _r[i])\n")
            write(f, "sweep = $(_sweep[i])\n")
            write(f, "dihed = $(_dihed[i])\n")
            write(f, "nelem_span = $(n_span)\n")
            write(f, "type_span = $(type_span)\n")
            write(f, "\n")
        end
    end
end
