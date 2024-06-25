function metal(T, n, p, sh, rh, ev, ph)
    # Props
    k_metal = 14.1091 + 0.0149964 * T   # [W/m-K]
    rho_metal = 8008.97 - 0.402716 * T  # [kg/m3]
    cp_metal = 501.933 + 0.196659 * T - 0.0000750219 * T^2  # [J/kg-K]
    
    # Discretize
    N_m1 = n.N_p
    N_m2a = N_m1 + n.N_sh
    N_m2b = N_m2a + n.N_sh
    N_m2c = N_m2b + n.N_ts
    N_m2d = N_m2c + n.N_ts
    N_m3a = N_m2d + n.N_rh
    N_m3b = N_m3a + n.N_rh
    N_m3c = N_m3b + n.N_ts
    N_m3d = N_m3c + n.N_ts
    N_m4a = N_m3d + n.N_ev
    N_m4b = N_m4a + n.N_ev
    N_m4c = N_m4b + n.N_ts
    N_m4d = N_m4c + n.N_ts
    N_m5a = N_m4d + n.N_ph
    N_m5b = N_m5a + n.N_ph
    N_m5c = N_m5b + n.N_ts
    N_m5d = N_m5c + n.N_ts
    
    # Pipe
    p.k_mw = k_metal[1:N_m1]
    p.rho_mw = rho_metal[1:N_m1]
    p.cp_mw = cp_metal[1:N_m1]
    
    # Superheater
    sh.k_m = k_metal[N_m1+1:N_m2a]
    sh.rho_m = rho_metal[N_m1+1:N_m2a]
    sh.cp_m = cp_metal[N_m1+1:N_m2a]
    sh.k_mw = k_metal[N_m2a+1:N_m2b]
    sh.rho_mw = rho_metal[N_m2a+1:N_m2b]
    sh.cp_mw = cp_metal[N_m2a+1:N_m2b]
    
    # ts
    sh.rho_uts = rho_metal[N_m2b+1:N_m2c]
    sh.k_uts = k_metal[N_m2b+1:N_m2c]
    sh.cp_uts = cp_metal[N_m2b+1:N_m2c]
    sh.rho_lts = rho_metal[N_m2c+1:N_m2d]
    sh.k_lts = k_metal[N_m2c+1:N_m2d]
    sh.cp_lts = cp_metal[N_m2c+1:N_m2d]
    
    # Reheater
    rh.k_m = k_metal[N_m2d+1:N_m3a]
    rh.rho_m = rho_metal[N_m2d+1:N_m3a]
    rh.cp_m = cp_metal[N_m2d+1:N_m3a]
    rh.k_mw = k_metal[N_m3a+1:N_m3b]
    rh.rho_mw = rho_metal[N_m3a+1:N_m3b]
    rh.cp_mw = cp_metal[N_m3a+1:N_m3b]
    
    # ts
    rh.rho_uts = rho_metal[N_m3b+1:N_m3c]
    rh.k_uts = k_metal[N_m3b+1:N_m3c]
    rh.cp_uts = cp_metal[N_m3b+1:N_m3c]
    rh.rho_lts = rho_metal[N_m3c+1:N_m3d]
    rh.k_lts = k_metal[N_m3c+1:N_m3d]
    rh.cp_lts = cp_metal[N_m3c+1:N_m3d]
    
    # Evaporator
    ev.k_m = k_metal[N_m3d+1:N_m4a]
    ev.rho_m = rho_metal[N_m3d+1:N_m4a]
    ev.cp_m = cp_metal[N_m3d+1:N_m4a]
    ev.k_mw = k_metal[N_m4a+1:N_m4b]
    ev.rho_mw = rho_metal[N_m4a+1:N_m4b]
    ev.cp_mw = cp_metal[N_m4a+1:N_m4b]
    
    # ts
    ev.rho_uts = rho_metal[N_m4b+1:N_m4c]
    ev.k_uts = k_metal[N_m4b+1:N_m4c]
    ev.cp_uts = cp_metal[N_m4b+1:N_m4c]
    ev.rho_lts = rho_metal[N_m4c+1:N_m4d]
    ev.k_lts = k_metal[N_m4c+1:N_m4d]
    ev.cp_lts = cp_metal[N_m4c+1:N_m4d]
    
    # Preheater
    ph.k_m = k_metal[N_m4d+1:N_m5a]
    ph.rho_m = rho_metal[N_m4d+1:N_m5a]
    ph.cp_m = cp_metal[N_m4d+1:N_m5a]
    ph.k_mw = k_metal[N_m5a+1:N_m5b]
    ph.rho_mw = rho_metal[N_m5a+1:N_m5b]
    ph.cp_mw = cp_metal[N_m5a+1:N_m5b]
    
    # ts
    ph.rho_uts = rho_metal[N_m5b+1:N_m5c]
    ph.k_uts = k_metal[N_m5b+1:N_m5c]
    ph.cp_uts = cp_metal[N_m5b+1:N_m5c]
    ph.rho_lts = rho_metal[N_m5c+1:N_m5d]
    ph.k_lts = k_metal[N_m5c+1:N_m5d]
    ph.cp_lts = cp_metal[N_m5c+1:N_m5d]
    
    return p, sh, rh, ev, ph
end
