# Import necessary packages
using SteamGen

dates = ["yyyy_mm_dd"]   # example day

# Simulation directives
@kwdef struct directives
    day::String = dates[1]
    start::Int64 = 1
    t_f::Int64 = 80800                                  # duration
    stop::Int64 = 80801                                 # start time + duration + 1
    step::Int64 = 80800                                 # number of steps
    tspan::Vector{Float64} = range(0, stop=80800, length=80800)   # time vector
    L::Float64 = 30.0                                       # length of p
    R_contact::Float64 = 0.0
    reduc::Float64 = 1.0                                    # reduction of heat transfer coefficient on tubesheets
    elbow::Float64 = 3.156                                # length of elbow
    temp::Float64 = 0.0                                     # run kelvin or celsius
    X::Float64 = 0.3                                      # vapor fraction (quality)
    T_offset::Float64 = 1.0
end
d = directives()

# Plot directives
pl = Dict(
    :plots => true,         # display any plots at all
    :save => false,         # save selected plots
    :bcs => true,           # display boundary conditions
    :temps => true          # display predicted temperature profiles
)

# Placeholder functions for geometry, crescentdunesdata, setup, solvedTdt, errors, plotstoplot
p, sh, rh, ev, ph, d = geometry(d);
sh, rh, ev, ph, d = assign_data(p, sh, rh, ev, ph, d)
n, p, sh, rh, ev, ph, T_i, d = setup(p, sh, rh, ev, ph, d)

[t,T] = ode15s(@(t,T) solvedTdt(T,t,n,p,sh,rh,ev,ph,d), d.tspan, T_i);
# T, n, p, sh, rh, ev, ph, d = errors(T, n, p, sh, rh, ev, ph, d)

# T, t, n, p, sh, rh, ev, ph, d, pl = plotstoplot(T, t, n, p, sh, rh, ev, ph, d, pl)

# # Crescent dunes data
# p, sh, rh, ev, ph, d = geometry(d)
# p[:L] = d[:L]
# sh[:R_contact] = d[:R_contact]
# rh[:R_contact] = d[:R_contact]
# ev[:R_contact] = d[:R_contact]
# ph[:R_contact] = d[:R_contact]

# (p, sh, rh, ev, ph, d) = crescentdunesdata(p, sh, rh, ev, ph, d)

# # ODE setup
# (n, p, sh, rh, ev, ph, ts, T_i, d) = setup(p, sh, rh, ev, ph, d)
# prob = ODEProblem((t, T) -> solvedTdt(T, t, n, p, sh, rh, ev, ph, d), T_i, d[:tspan])
# sol = solve(prob, Rodas5())

# # Error
# (T, n, p, sh, rh, ev, ph, d) = errors(sol.u, n, p, sh, rh, ev, ph, d)

# # PLOTS
# if pl[:plots]
#     (T, t, n, p, sh, rh, ev, ph, d, pl) = plotstoplot(T, sol.t, n, p, sh, rh, ev, ph, d, pl)
# else
#     println("SPEED TEST BABY!!! BRRRRRRRR")
# end
