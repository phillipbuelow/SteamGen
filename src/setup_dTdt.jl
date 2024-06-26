
function solvedTdt(T, t, n, p, sh, rh, ev, ph, d)

    dTdt = zeros(n.N_lts_ph,1)       # initialize dTdt vector
    
    ## Extract temperatures from T vector
    T_s_p    = T[1:n.N_s_p]                  # Pipe temperatures in primary heat exchanger
    T_mw_p   = T[n.N_s_p+1:n.N_mw_p]         # Metal wall temperatures in primary heat exchanger
    T_s_sh   = T[n.N_mw_p+1:n.N_s_sh]        # Salt temperatures in superheater
    T_m_sh   = T[n.N_s_sh+1:n.N_m_sh]        # Metal wall temperatures in superheater
    T_w_sh   = T[n.N_m_sh+1:n.N_w_sh]        # Water temperatures in superheater
    T_mw_sh  = T[n.N_w_sh+1:n.N_mw_sh]       # Metal wall temperatures in superheater
    T_uts_sh = T[n.N_mw_sh+1:n.N_uts_sh]     # Upper tubesheet temperatures in superheater
    T_lts_sh = T[n.N_uts_sh+1:n.N_lts_sh]    # Lower tubesheet temperatures in superheater
    T_s_rh   = T[n.N_lts_sh+1:n.N_s_rh]       # Salt temperatures in reheater
    T_m_rh   = T[n.N_s_rh+1:n.N_m_rh]         # Metal wall temperatures in reheater
    T_w_rh   = T[n.N_m_rh+1:n.N_w_rh]         # Water temperatures in reheater
    T_mw_rh  = T[n.N_w_rh+1:n.N_mw_rh]        # Metal wall temperatures in reheater
    T_uts_rh = T[n.N_mw_rh+1:n.N_uts_rh]      # Upper tubesheet temperatures in reheater
    T_lts_rh = T[n.N_uts_rh+1:n.N_lts_rh]     # Lower tubesheet temperatures in reheater
    T_s_ev   = T[n.N_lts_rh+1:n.N_s_ev]       # Salt temperatures in evaporator
    T_m_ev   = T[n.N_s_ev+1:n.N_m_ev]         # Metal wall temperatures in evaporator
    T_w_ev   = T[n.N_m_ev+1:n.N_w_ev]         # Water temperatures in evaporator
    T_mw_ev  = T[n.N_w_ev+1:n.N_mw_ev]        # Metal wall temperatures in evaporator
    T_uts_ev = T[n.N_mw_ev+1:n.N_uts_ev]      # Upper tubesheet temperatures in evaporator
    T_lts_ev = T[n.N_uts_ev+1:n.N_lts_ev]     # Lower tubesheet temperatures in evaporator
    T_s_ph   = T[n.N_lts_ev+1:n.N_s_ph]       # Salt temperatures in preheater
    T_m_ph   = T[n.N_s_ph+1:n.N_m_ph]         # Metal wall temperatures in preheater
    T_w_ph   = T[n.N_m_ph+1:n.N_w_ph]         # Water temperatures in preheater
    T_mw_ph  = T[n.N_w_ph+1:n.N_mw_ph]        # Metal wall temperatures in preheater
    T_uts_ph = T[n.N_mw_ph+1:n.N_uts_ph]      # Upper tubesheet temperatures in preheater
    T_lts_ph = T[n.N_uts_ph+1:n.N_lts_ph]     # Lower tubesheet temperatures in preheater
    
    ## Calculate dTdt for each component
    
    # Primary heat exchanger (pipe]
    dTdt[1] = (p.h_s[1] * p.A_s * (T_mw_p[1] - T_s_p[1]) + sh.mdot_s[tt] * bcs.cp_s[1] * sh.T_s_in[tt] - sh.mdot_s[tt] * p.cp_s[1] * T_s_p[1]) / (p.rho_s[1] * p.cp_s[1] * p.V_s)
    dTdt[2:n.N_s_p] = (p.h_s[2:end] .* p.A_s .* (T_mw_p[2:end] - T_s_p[2:end]) + sh.mdot_s[tt] 
                        .* p.cp_s[1:end-1].* T_s_p[1:end-1] - sh.mdot_s[tt] .* p.cp_s[2:end] 
                        .* T_s_p[2:end]) ./ (p.rho_s[2:end] .* p.cp_s[2:end] .* p.V_s)
    
    # Primary heat exchanger (metal wall)
    dTdt[n.N_s_p+1] = ((T_s_p[1] - T_mw_p[1]) ./ (1 ./ (p.h_s[1] * p.A_s) + log(p.r_c_p / p.r_i_p) / (2 * pi * p.dx * p.k_mw[1]))
        + (sh.T_amb - T_mw_p[1]) ./ (1 ./ (sh.h_amb * p.A_s_insul) + log(p.r_o_insul / p.r_o_p) / (2 * pi * p.dx * k_insul[1]) + log(p.r_o_p / p.r_c_p) / (2 * pi * p.dx * p.k_mw[1]))
        + (T_mw_p[2] - T_mw_p[1]) * p.k_mw[1] * p.A_c_mw / p.dx) ./ (p.rho_mw[1] * p.cp_mw[1] * p.V_mw)
    dTdt[n.N_s_p+2:n.N_mw_p-1] = ((T_s_p[2:end-1] - T_mw_p[2:end-1]) ./ (1 ./ (p.h_s[2:end-1] * p.A_s) + log(p.r_c_p / p.r_i_p) / (2 * pi * p.dx * p.k_mw[2:end-1])) 
        + (sh.T_amb - T_mw_p[2:end-1]) ./ (1 ./ (sh.h_amb * p.A_s_insul) + log(p.r_o_insul / p.r_o_p) / (2 * pi * p.dx * k_insul[2:end-1]) + log(p.r_o_p / p.r_c_p) / (2 * pi * p.dx * p.k_mw[2:end-1])) 
        + (T_mw_p[1:end-2] - T_mw_p[2:end-1]) * p.k_mw[2:end-1] * p.A_c_mw / p.dx 
        + (T_mw_p[3:end] - T_mw_p[2:end-1]) * p.k_mw[2:end-1] * p.A_c_mw / p.dx) ./ (p.rho_mw[2:end-1] * p.cp_mw[2:end-1] * p.V_mw)

    dTdt[n.N_mw_p] = ((T_s_p[end] - T_mw_p[end]) ./ (1 ./ (p.h_s[end] * p.A_s) + log(p.r_c_p / p.r_i_p) / (2 * pi * p.dx * p.k_mw[end])) 
        + (sh.T_amb - T_mw_p[end]) ./ (1 ./ (sh.h_amb * p.A_s_insul) + log(p.r_o_insul / p.r_o_p) / (2 * pi * p.dx * k_insul[end]) + log(p.r_o_p / p.r_c_p) / (2 * pi * p.dx * p.k_mw[end])) 
        + (T_mw_p[end-1] - T_mw_p[end]) * p.k_mw[end] * p.A_c_mw / p.dx) ./ (p.rho_mw[end] * p.cp_mw[end] * p.V_mw)
    
    # Superheater (salt)
    dTdt[n.N_mw_p+1] = (sh.h_s[1] * sh.A_s * (T_m_sh[1] - T_s_sh[1]) + sh.mdot_s[tt] * p.cp_s[end] * T_s_p[end] - sh.mdot_s[tt] * sh.cp_s[1] * T_s_sh[1]) / (sh.rho_s[1] * sh.cp_s[1] * sh.V_s)
    dTdt[n.N_mw_p+2:n.N_s_sh] = (sh.h_s[2:end] * sh.A_s .* (T_m_sh[2:end] - T_s_sh[2:end]) + sh.h_s[2:end] * sh.A_mw .* (T_mw_sh[2:end] - T_s_sh[2:end]) 
        + sh.mdot_s[tt] * sh.cp_s[1:end-1] .* T_s_sh[1:end-1] - sh.mdot_s[tt] * sh.cp_s[2:end] .* T_s_sh[2:end]) ./ (sh.rho_s[2:end] * sh.cp_s[2:end] * sh.V_s)
    
    # Superheater (metal wall)
    dTdt[n.N_s_sh+1] = ((T_s_sh[1] - T_m_sh[1]) ./ (1 ./ (sh.h_s[1] * sh.A_s) + log(sh.r_c_sh / sh.r_i_sh) / (2 * pi * sh.dx * sh.k_mw[1])) 
        + (sh.T_amb - T_m_sh[1]) ./ (1 ./ (sh.h_amb * sh.A_s_insul) + log(sh.r_o_insul / sh.r_o_sh) / (2 * pi * sh.dx * k_insul_sh[1]) + log(sh.r_o_sh / sh.r_c_sh) / (2 * pi * sh.dx * sh.k_mw[1])) 
        + (T_m_sh[2] - T_m_sh[1]) * sh.k_mw[1] * sh.A_c_mw / sh.dx) ./ (sh.rho_mw[1] * sh.cp_mw[1] * sh.V_mw)
    dTdt[n.N_s_sh+2:n.N_m_sh-1] = ((T_s_sh[2:end-1] - T_m_sh[2:end-1]) ./ (1 ./ (sh.h_s[2:end-1] * sh.A_s) + log(sh.r_c_sh / sh.r_i_sh) / (2 * pi * sh.dx * sh.k_mw[2:end-1])) 
        + (sh.T_amb - T_m_sh[2:end-1]) ./ (1 ./ (sh.h_amb * sh.A_s_insul) + log(sh.r_o_insul / sh.r_o_sh) / (2 * pi * sh.dx * k_insul_sh[2:end-1]) + log(sh.r_o_sh / sh.r_c_sh) / (2 * pi * sh.dx * sh.k_mw[2:end-1])) 
        + (T_m_sh[1:end-2] - T_m_sh[2:end-1]) * sh.k_mw[2:end-1] * sh.A_c_mw / sh.dx 
        + (T_m_sh[3:end] - T_m_sh[2:end-1]) * sh.k_mw[2:end-1] * sh.A_c_mw / sh.dx) ./ (sh.rho_mw[2:end-1] * sh.cp_mw[2:end-1] * sh.V_mw)
    dTdt[n.N_m_sh] = ((T_s_sh[end] - T_m_sh[end]) ./ (1 ./ (sh.h_s[end] * sh.A_s) + log(sh.r_c_sh / sh.r_i_sh) / (2 * pi * sh.dx * sh.k_mw[end])) 
        + (sh.T_amb - T_m_sh[end]) ./ (1 ./ (sh.h_amb * sh.A_s_insul) + log(sh.r_o_insul / sh.r_o_sh) / (2 * pi * sh.dx * k_insul_sh[end]) + log(sh.r_o_sh / sh.r_c_sh) / (2 * pi * sh.dx * sh.k_mw[end])) 
        + (T_m_sh[end-1] - T_m_sh[end]) * sh.k_mw[end] * sh.A_c_mw / sh.dx) ./ (sh.rho_mw[end] * sh.cp_mw[end] * sh.V_mw)
    
    # Superheater (water)
    dTdt[n.N_m_sh+1] = (sh.mdot_w[tt] * bcs.cp_w[tt] * sh.T_w_in[tt] - sh.mdot_w[tt] * sh.cp_w[1] * T_w_sh[1]) / (sh.rho_w[1] * sh.cp_w[1] * sh.V_w)
    dTdt[n.N_m_sh+2:n.N_w_sh] = (sh.mdot_w[tt] * sh.cp_w[1:end-1] .* T_w_sh[1:end-1] - sh.mdot_w[tt] * sh.cp_w[2:end] .* T_w_sh[2:end]) ./ (sh.rho_w[2:end] * sh.cp_w[2:end] * sh.V_w)
    
    return dTdt
end