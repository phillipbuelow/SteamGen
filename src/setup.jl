function setup(p, sh, rh, ev, ph, d)
    # Initialize n struct
    n = Dict()

    # Discretization
    n[:N_p] = 50
    n[:N_s_p] = Int64(n[:N_p])
    n[:N_mw_p] = Int64(n[:N_s_p] + n[:N_p])

    # Superheater
    n[:N_sh] = Int64(sh[:N_baffles])
    n[:N_s_sh] = Int64(n[:N_mw_p] + n[:N_sh])
    n[:N_m_sh] = Int64(n[:N_s_sh] + n[:N_sh])
    n[:N_w_sh] = Int64(n[:N_m_sh] + n[:N_sh])
    n[:N_mw_sh] = Int64(n[:N_w_sh] + n[:N_sh])
    n[:N_ts] = 20
    n[:N_uts_sh] = Int64(n[:N_mw_sh] + n[:N_ts])
    n[:N_lts_sh] = Int64(n[:N_uts_sh] + n[:N_ts])

    # Reheater
    n[:N_rh] = Int64(rh[:N_baffles])
    n[:N_s_rh] = Int64(n[:N_lts_sh] + n[:N_rh])
    n[:N_m_rh] = Int64(n[:N_s_rh] + n[:N_rh])
    n[:N_w_rh] = Int64(n[:N_m_rh] + n[:N_rh])
    n[:N_mw_rh] = Int64(n[:N_w_rh] + n[:N_rh])
    n[:N_uts_rh] = Int64(n[:N_mw_rh] + n[:N_ts])
    n[:N_lts_rh] = Int64(n[:N_uts_rh] + n[:N_ts])

    # Evaporator
    n[:N_ev] = Int64(ev[:N_baffles])
    n[:N_s_ev] = Int64(n[:N_lts_rh] + n[:N_ev])
    n[:N_m_ev] = Int64(n[:N_s_ev] + n[:N_ev])
    n[:N_w_ev] = Int64(n[:N_m_ev] + n[:N_ev])
    n[:N_mw_ev] = Int64(n[:N_w_ev] + n[:N_ev])
    n[:N_uts_ev] = Int64(n[:N_mw_ev] + n[:N_ts])
    n[:N_lts_ev] = Int64(n[:N_uts_ev] + n[:N_ts])

    # Preheater
    n[:N_ph] = Int64(ph[:N_baffles])
    n[:N_s_ph] = Int64(n[:N_lts_ev] + n[:N_ph])
    n[:N_m_ph] = Int64(n[:N_s_ph] + n[:N_ph])
    n[:N_w_ph] = Int64(n[:N_m_ph] + n[:N_ph])
    n[:N_mw_ph] = Int64(n[:N_w_ph] + n[:N_ph])
    n[:N_uts_ph] = Int64(n[:N_mw_ph] + n[:N_ts])
    n[:N_lts_ph] = Int64(n[:N_uts_ph] + n[:N_ts])

    # Discretize x-direction
    @show n[:N_sh]
    p[:x] = range(0, stop=p[:L], length=n[:N_p])
    sh[:x] = range(0, stop=sh[:L], length=n[:N_sh])
    rh[:x] = range(0, stop=rh[:L], length=n[:N_rh])
    ev[:x] = range(0, stop=ev[:L], length=n[:N_ev])
    ph[:x] = range(0, stop=ph[:L], length=n[:N_ph])

    # Initial temperature
    T_i = zeros(Float64, sum(values(n)))

    # Set initial temperatures for each section
    idx = 1

    # p
    @show n[:N_s_p]
    @show idx
    T_i[idx:(idx+n[:N_s_p]-1)] .= sh[:T_s_in][1]
    idx += Int64(n[:N_s_p])
    T_i[idx:idx+n[:N_mw_p]-1] .= sh[:T_s_in][1]
    idx += Int64(n[:N_mw_p])

    # superheater
    T_i[idx:idx+n[:N_sh]-1] = range(sh[:T_s_in][1], sh[:T_s_out][1], length=n[:N_sh])
    idx += Int64(n[:N_sh])
    T_i[idx:idx+n[:N_sh]-1] = range(sh[:T_w_out][1], sh[:T_w_in][1], length=n[:N_sh])
    idx += Int64(n[:N_sh])
    T_i[idx:idx+n[:N_sh]-1] = (T_i[idx-n[:N_sh]:idx-1] .+ T_i[idx-n[:N_sh]:idx-1]) ./ 2
    idx += Int64(n[:N_sh])
    T_i[idx:idx+n[:N_sh]-1] = T_i[idx-n[:N_sh]:idx-1]
    idx += Int64(n[:N_sh])
    T_i[idx:idx+n[:N_ts]-1] = range(sh[:T_s_in][1], sh[:T_w_out][1], length=n[:N_ts])
    idx += Int64(n[:N_ts])
    T_i[idx:idx+n[:N_ts]-1] = fill(sh[:T_w_in][1], n[:N_ts])
    idx += Int64(n[:N_ts])

    # reheater
    T_i[idx:idx+n[:N_rh]-1] = range(rh[:T_s_in][1], rh[:T_s_out][1], length=n[:N_rh])
    idx += Int64(n[:N_rh])
    T_i[idx:idx+n[:N_rh]-1] = range(sh[:T_w_out][1], rh[:T_w_in][1], length=n[:N_rh])
    idx += Int64(n[:N_rh])
    T_i[idx:idx+n[:N_rh]-1] = (T_i[idx-n[:N_rh]:idx-1] .+ T_i[idx-n[:N_rh]:idx-1]) ./ 2
    idx += Int64(n[:N_rh])
    T_i[idx:idx+n[:N_rh]-1] = T_i[idx-n[:N_rh]:idx-1]
    idx += Int64(n[:N_rh])
    T_i[idx:idx+n[:N_ts]-1] = range(rh[:T_s_in][1], sh[:T_w_out][1], length=n[:N_ts])
    idx += Int64(n[:N_ts])
    T_i[idx:idx+n[:N_ts]-1] = fill(rh[:T_w_in][1], n[:N_ts])
    idx += Int64(n[:N_ts])

    # evaporator
    T_i[idx:idx+n[:N_ev]-1] = range(ev[:T_s_in][1], ev[:T_s_out][1], length=n[:N_ev])
    idx += Int64(n[:N_ev])
    T_i[idx:idx+n[:N_ev]-1] = range(ev[:T_w_out][1], ev[:T_w_in][1], length=n[:N_ev])
    idx += Int64(n[:N_ev])
    T_i[idx:idx+n[:N_ev]-1] = (T_i[idx-n[:N_ev]:idx-1] .+ T_i[idx-n[:N_ev]:idx-1]) ./ 2
    idx += Int64(n[:N_ev])
    T_i[idx:idx+n[:N_ev]-1] = T_i[idx-n[:N_ev]:idx-1]
    idx += Int64(n[:N_ev])
    T_i[idx:idx+n[:N_ts]-1] = range(ev[:T_s_in][1], ev[:T_w_out][1], length=n[:N_ts])
    idx += Int64(n[:N_ts])
    T_i[idx:idx+n[:N_ts]-1] = fill(ev[:T_w_in][1], n[:N_ts])
    idx += Int64(n[:N_ts])

    # preheater
    T_i[idx:idx+n[:N_ph]-1] = range(ph[:T_s_in][1], ph[:T_s_out][1], length=n[:N_ph])
    idx += Int64(n[:N_ph])
    T_i[idx:idx+n[:N_ph]-1] = range(ph[:T_w_out][1], ph[:T_w_in][1], length=n[:N_ph])
    idx += Int64(n[:N_ph])
    T_i[idx:idx+n[:N_ph]-1] = (T_i[idx-n[:N_ph]:idx-1] .+ T_i[idx-n[:N_ph]:idx-1]) ./ 2
    idx += Int64(n[:N_ph])
    T_i[idx:idx+n[:N_ph]-1] = T_i[idx-n[:N_ph]:idx-1]
    idx += Int64(n[:N_ph])
    T_i[idx:idx+n[:N_ts]-1] = range(ph[:T_s_in][1], ph[:T_w_out][1], length=n[:N_ts])
    idx += Int64(n[:N_ts])
    T_i[idx:idx+n[:N_ts]-1] = fill(ph[:T_w_in][1], n[:N_ts])
    idx += Int64(n[:N_ts])

    # Pressure matrices
    sh[:P] = ones(length(sh[:P_w_in]), n[:N_sh])
    rh[:P] = ones(length(rh[:P_w_in]), n[:N_rh])
    ev[:P] = ones(length(ev[:P_w_in]), n[:N_ev])
    ph[:P] = ones(length(ph[:P_w_in]), n[:N_ph])

    # Discretize p
    sh[:p] = range(0, stop=sh[:P_w_out][1], length=n[:N_sh])
    rh[:p] = range(0, stop=rh[:P_w_out][1], length=n[:N_rh])
    ev[:p] = range(0, stop=ev[:P_w_out][1], length=n[:N_ev])
    ph[:p] = range(0, stop=ph[:P_w_out][1], length=n[:N_ph])

    return n, p, sh, rh, ev, ph, T_i, d
end
