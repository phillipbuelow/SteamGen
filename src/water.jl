function water(bcw, sh, rh, ev, ph, T, n, tt, d)
    # Concatenate pressures
    P = vcat(sh.P[tt,:], rh.P[tt,:], ev.P[tt,:], ph.P[tt,:])

    # Assign bcw.P
    bcw.P = [sh.P_w_in[tt], rh.P_w_in[tt], ev.P_w_in[tt], ph.P_w_in[tt]]

    # Adjust temperatures
    bcw.T = [sh.T_w_in[tt], rh.T_w_in[tt], ev.T_w_in[tt], ph.T_w_in[tt]] - d.temp

    # Discretize sections
    N_f2 = n["N_sh"]
    N_f3 = N_f2 + n["N_rh"]
    N_f4 = N_f3 + n["N_ev"]
    N_f5 = N_f4 + n["N_ph"]

    # Properties calculation
    bcw.cp_w = similar(bcw.T)
    bcw.cp_w[1] = XSteamI("cp_pT", bcw.P[1] - 0.96, bcw.T[1]) * 1000
    bcw.cp_w[2] = XSteamI("cp_pT", bcw.P[2], bcw.T[2]) * 1000
    bcw.cp_w[3] = XSteamI("cp_pT", bcw.P[3], bcw.T[3]) * 1000
    bcw.cp_w[4] = XSteamI("cp_pT", bcw.P[4], bcw.T[4]) * 1000

    # Evaporator properties (multiphase)
    cp_w_V = XSteamI("cpV_T", T[N_f3+1:N_f4])' .* 1000
    k_w_V = XSteamI("tcV_T", T[N_f3+1:N_f4])'
    rho_w_V = XSteamI("rhoV_T", T[N_f3+1:N_f4])'

    # Evaporator properties (single phase)
    cp_w_L = XSteamI("cpL_T", T[N_f3+1:N_f4])' .* 1000
    k_w_L = XSteamI("tcL_T", T[N_f3+1:N_f4])'
    rho_w_L = XSteamI("rhoL_T", T[N_f3+1:N_f4])'

    # Assign properties to evaporator
    ev.cp_w = (1.0 - ev.X) .* cp_w_L + ev.X .* cp_w_V
    ev.k_w = (1.0 - ev.X) .* k_w_L + ev.X .* k_w_V
    ev.rho_w = (1.0 - ev.X) .* rho_w_L + ev.X .* rho_w_V
    ev.mu_w = 0.00002  # Example value, adjust as necessary

    # Superheater properties
    sh.cp_w = XSteamI("cp_pT", P[1:N_f2]', T[1:N_f2]')' .* 1000
    sh.k_w = XSteamI("tc_pT", P[1:N_f2]', T[1:N_f2]')'
    sh.rho_w = XSteamI("rho_pT", P[1:N_f2]', T[1:N_f2]')'
    sh.mu_w = XSteamI("my_pT", P[1:N_f2]', T[1:N_f2]')'
    sh.Pr_w = (sh.mu_w .* sh.cp_w) ./ sh.k_w
    sh.nu_w = sh.mu_w ./ sh.rho_w

    # Reheater properties
    rh.cp_w = XSteamI("cp_pT", P[N_f2+1:N_f3]', T[N_f2+1:N_f3]')' .* 1000
    rh.k_w = XSteamI("tc_pT", P[N_f2+1:N_f3]', T[N_f2+1:N_f3]')'
    rh.rho_w = XSteamI("rho_pT", P[N_f2+1:N_f3]', T[N_f2+1:N_f3]')'
    rh.mu_w = XSteamI("my_pT", P[N_f2+1:N_f3]', T[N_f2+1:N_f3]')'
    rh.Pr_w = (rh.mu_w .* rh.cp_w) ./ rh.k_w
    rh.nu_w = rh.mu_w ./ rh.rho_w

    # Preheater properties
    ph.cp_w = XSteamI("cp_pT", P[N_f4+1:N_f5]', T[N_f4+1:N_f5]')' .* 1000
    ph.k_w = XSteamI("tc_pT", P[N_f4+1:N_f5]', T[N_f4+1:N_f5]')'
    ph.rho_w = XSteamI("rho_pT", P[N_f4+1:N_f5]', T[N_f4+1:N_f5]')'
    ph.mu_w = XSteamI("my_pT", P[N_f4+1:N_f5]', T[N_f4+1:N_f5]')'
    ph.Pr_w = (ph.mu_w .* ph.cp_w) ./ ph.k_w
    ph.nu_w = ph.mu_w ./ ph.rho_w

    # Calculate convective coefficients
    calculate_convective_coefficients(sh, tt, n)
    calculate_convective_coefficients(rh, tt, n)
    calculate_convective_coefficients(ev, tt, n)
    calculate_convective_coefficients(ph, tt, n)

    return bcw, sh, rh, ev, ph
end

function calculate_convective_coefficients(section, tt, n)
    u = (section.mdot_w[tt] / section.N_tubes) ./ (section.rho_w .* section.A_c_tube)
    L_c = section.D_i_t

    section.Re_w = u .* L_c ./ section.nu_w

    section.Nu = similar(section.Re_w)
    for i in 1:length(section.Re_w)
        if section.Re_w[i] < 2300
            section.Nu[i] = 1.86 * (section.Re_w[i] .* section.Pr_w[i] .* L_c ./ section.dx).^(1/3)
        elseif section.Re_w[i] >= 2300 && section.Re_w[i] <= 10000
            f = (100 * section.Re_w[i])^(-1/4)
            section.Nu[i] = (f / 8 * (section.Re_w[i] - 1000) * section.Pr_w[i]) / (1 + 12.7 * (f / 8)^(1/2) * (section.Pr_w[i]^(2/3) - 1))
        elseif section.Re_w[i] > 10000
            f = (1.58 * log(section.Re_w[i]) - 3.28)^(-2)
            section.Nu[i] = (f / 2 * section.Re_w[i] * section.Pr_w[i]) / (1.07 + 12.7 * (f / 2)^0.5 * (section.Pr_w[i]^(2/3) - 1))
        end
    end

    section.h_w = section.Nu .* section.k_w ./ L_c
end
