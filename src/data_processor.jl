function assign_data(p, sh, rh, ev, ph, d)
    directory = "C:/Users/phill/.julia/dev/SteamGen/data/transient_bcs/"
    data_file = joinpath(directory, string(d.day) * ".jld2")
    data = load(data_file)["data"]

    # superheater
    sh[:T_w_in]  = round.(data["stmfrmdrum1tosheattmp"][d.start:end] .+ d.T_offset, digits=1)
    sh[:T_w_out] = round.(data["sheat1outtemp"][d.start:end] .+ d.T_offset, digits=1)
    sh[:T_s_in]  = round.(data["mixedsalttosgs1tmpseltd"][d.start:end] .+ d.T_offset, digits=1)
    sh[:T_s_out] = round.(data["shtr1saltouttemp"][d.start:end] .+ d.T_offset, digits=1)
    sh[:mdot_w]  = round.(data["shtr1stm_stflowseltd"][d.start:end], digits=2)
    sh[:mdot_s]  = round.(data["shtr1saltoutflow"][d.start:end], digits=2)
    sh[:P_w_in]  = round.(data["stmfrmdrum1tosheatpress"][d.start:end] .* 0.01, digits=2)
    sh[:P_w_out] = round.(data["shtr1stm_stpress"][d.start:end] .* 0.01, digits=2)

    for i in 1:d.stop
        if sh[:mdot_w][i] < 1
            sh[:mdot_w][i] = 1
        end
    end

    dlength = d.start + length(sh[:T_w_in]) - 1

    # reheater
    rh[:T_s_in]   = round.(data["mixedsalttosgs1tmpseltd"][d.start:end] .+ d.T_offset, digits=1)
    rh[:T_s_out]  = round.(data["rhtr1saltouttmp"][d.start:end] .+ d.T_offset, digits=1)
    rh[:mdot_s]   = round.(data["rhtr1saltoutflow"][d.start:end], digits=2)
    rh[:mdot_w]   = fill(0.0, dlength)
    rh[:T_w_in]   = fill(0.0, dlength)
    rh[:P_w_in]   = fill(0.0, dlength)
    rh[:P_w_drop] = fill(0.0, dlength)

    for i in 1:length(rh[:T_s_in])
        if sh[:mdot_w][i] < 6.5
            rh[:mdot_w][i]   = sh[:mdot_w][i]
            rh[:P_w_in][i]   = sh[:P_w_out][i] * 0.07
            rh[:P_w_drop][i] = (sh[:P_w_in][i] - sh[:P_w_out][i]) * 0.07
            rh[:T_w_in][i]   = sh[:T_w_out][i]
        elseif 6.5 <= sh[:mdot_w][i] && sh[:mdot_w][i] <= 30
            rh[:mdot_w][i]   = sh[:mdot_w][i] * 0.75
            rh[:P_w_in][i]   = sh[:P_w_out][i] * 0.15
            rh[:P_w_drop][i] = (sh[:P_w_in][i] - sh[:P_w_out][i]) * 0.15
            rh[:T_w_in][i]   = sh[:T_w_out][i] * 0.8
        elseif sh[:mdot_w][i] > 30
            rh[:mdot_w][i]   = sh[:mdot_w][i] * 0.8
            rh[:P_w_in][i]   = sh[:P_w_out][i] * 0.2
            rh[:P_w_drop][i] = (sh[:P_w_in][i] - sh[:P_w_out][i]) * 0.2
            rh[:T_w_in][i]   = sh[:T_w_out][i] * 0.6
        end
    end

    rh[:P_w_out] = rh[:P_w_in] .- rh[:P_w_drop]

    # evaporator
    ev[:T_w_in]  = round.(data["evapor1todrum1tempseltd"][d.start:end], digits=1)
    ev[:T_w_out] = round.(data["evapor1todrum1tempseltd"][d.start:end], digits=1)
    ev[:P_w_in]  = round.(data["ev_recircpmps_press"][d.start:end] .* 0.01, digits=2)
    ev[:P_w_out] = round.(data["ev_todrum_press"][d.start:end] .* 0.01, digits=2)
    ev[:mdot_w]  = round.(data["sgsrecircpmps_mdot"][d.start:end], digits=1)
    ev[:mdot_s]  = sh[:mdot_s] .+ rh[:mdot_s]
    ev[:T_s_in]  = round.(data["mixedsalttoevap1tempsel"][d.start:end], digits=1)
    ev[:T_s_out] = round.(data["evap1saltotltmpseltd"][d.start:end], digits=1)
    ev[:X]       = d.X

    for i in 1:d.stop
        if ev[:mdot_w][i] < 1
            ev[:mdot_w][i] = 1
        end
    end

    # preheater
    ph[:T_s_in]  = round.(data["evap1saltotltmpseltd"][d.start:end], digits=1)
    ph[:T_s_out] = round.(data["econ1saltotltmpseltd"][d.start:end], digits=1)
    ph[:mdot_s]  = sh[:mdot_s] .+ rh[:mdot_s]
    ph[:mdot_w]  = sh[:mdot_w]
    ph[:P_w_in]  = 145 .* ones(dlength)
    ph[:P_w_out] = 145 .* ones(dlength)
    ph[:T_w_in] = ev[:T_w_out] * 0.80

    for i in 1:dlength
        if ph[:mdot_w][i] < 1
            ph[:mdot_w][i] = 1
        end
        if sh[:mdot_w][i] > 30
            ph[:T_w_in][i] = ev[:T_w_out][i] * 0.93
        elseif 6.5 < sh[:mdot_w][i] <= 30
            ph[:T_w_in][i] = ev[:T_w_out][i] * 0.90
        elseif sh[:mdot_w][i] <= 6.5
            ph[:T_w_in] = ev[:T_w_out] * 0.80
        end
    end
    ph[:T_w_out] = ev[:T_w_in]  # ASSUMPTION

    # ambient
    sh[:T_amb] = 25
    sh[:h_amb] = 0  # [W/m2-K]

    return sh, rh, ev, ph, d
end
