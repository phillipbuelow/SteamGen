# Import necessary packages
using Dates
using LinearAlgebra
using DifferentialEquations
using Plots

# Simulation directives
dates = ["01-10"] # valid days
d = Dict(
    :day => dates[1],
    :start => 1,
    :t_f => 80800, # duration
    :stop => 80801, # start time + duration + 1
    :step => 80800, # number of steps
    :tspan => range(0, stop=80800, length=80800), # time vector
    :L => 30, # length of p
    :R_contact => 0,
    :reduc => 1, # reduction of heat transfer coefficient on tubesheets
    :elbow => 3.156, # length of elbow
    :temp => 0, # run kelvin or celsius
    :X => 0.3 # vapor fraction (quality)
)

# Plot directives
pl = Dict(
    :plots => true, # display any plots at all
    :save => false, # save selected plots
    :bcs => true, # display boundary conditions
    :temps => true # display predicted temperature profiles
)

# Placeholder functions for geometry, crescentdunesdata, setup, solvedTdt, errors, plotstoplot
function geometry(d)
    return (p, sh, rh, ev, ph, d)
end

function crescentdunesdata(p, sh, rh, ev, ph, d)
    return (p, sh, rh, ev, ph, d)
end

function setup(p, sh, rh, ev, ph, d)
    return (n, p, sh, rh, ev, ph, ts, T_i, d)
end

function solvedTdt(T, t, n, p, sh, rh, ev, ph, d)
    # Define your ODE system here
end

function errors(T, n, p, sh, rh, ev, ph, d)
    return (T, n, p, sh, rh, ev, ph, d)
end

function plotstoplot(T, t, n, p, sh, rh, ev, ph, d, pl)
    return (T, t, n, p, sh, rh, ev, ph, d, pl)
end

# Crescent dunes data
(p, sh, rh, ev, ph, d) = geometry(d)
p[:L] = d[:L]
sh[:R_contact] = d[:R_contact]
rh[:R_contact] = d[:R_contact]
ev[:R_contact] = d[:R_contact]
ph[:R_contact] = d[:R_contact]

(p, sh, rh, ev, ph, d) = crescentdunesdata(p, sh, rh, ev, ph, d)

# ODE setup
(n, p, sh, rh, ev, ph, ts, T_i, d) = setup(p, sh, rh, ev, ph, d)
prob = ODEProblem((t, T) -> solvedTdt(T, t, n, p, sh, rh, ev, ph, d), T_i, d[:tspan])
sol = solve(prob, Rodas5())

# Error
(T, n, p, sh, rh, ev, ph, d) = errors(sol.u, n, p, sh, rh, ev, ph, d)

# PLOTS
if pl[:plots]
    (T, t, n, p, sh, rh, ev, ph, d, pl) = plotstoplot(T, sol.t, n, p, sh, rh, ev, ph, d, pl)
else
    println("SPEED TEST BABY!!! BRRRRRRRR")
end
