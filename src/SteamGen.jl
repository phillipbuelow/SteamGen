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

## Generates globally used functions
include("structs.jl")
include("geometry.jl")
include("data_processor.jl")
include("setup.jl")
include("setup_dTdt.jl")

# Exportation of global functions 
export geometry
export assign_data
export setup
export solvedTdt

end
