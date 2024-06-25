module SteamGen

using Revise
using Dates
using Interpolations
using JLD2
using DataFrames
using LinearAlgebra
using DifferentialEquations
using Plots
using MAT

import Base: @kwdef

## Systems, general
include("structs.jl")
include("geometry.jl")

end
