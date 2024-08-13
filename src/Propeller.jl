"""
Propeller model based on blade element momentum theory using a linear inflow model.
"""
module Propeller

using FLOWMath
using PyFormattedStrings
using Dierckx
import AirfoilPolars as AP
using Random
using GridSpacing
using LinearAlgebra
using Statistics
using Trapz
using PyFormattedStrings
using Base.Threads
using ProgressMeter
using PrettySections
import Optim
using Interpolations
using FlightConditions
import GeometricTools as gt
using NumericalIntegration
import AirfoilFast as AF
using DataFrames
using CSV


include("polars.jl")
export PolarGrid, polar_interpolators, copy, correction3D
# export plot


include("hubtip_correction.jl")
export hubtip_correction_custom, hubtip_correction_none, hubtip_correction_modprandtl, hubtip_correction_prandtl


include("rotor.jl")
export span2radius, radius2span, radial_segments, span, r_hub, midpoints, linear_subdivide
export Rotor, correction3D!
export plot 
export update_chord!, update_twist!
export copy


include("bem.jl")
export rpm2thrust
export simulate_rpm2thrust, AeroResult


include("io.jl")
export export_flowunsteady


include("grid.jl")
export blade2grid, map_loads

include("estimations.jl")
export estimate_Re

end
