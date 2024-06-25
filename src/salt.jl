function salt(T, p, sh, rh, ev, ph, tt, n, d)
    bcs = Dict{String, Any}()
    
    # Convert structs to dictionaries
    p = Dict{Symbol, Any}(Symbol(name) => Dict{String, Any}() for name in fieldnames(p))
    sh = Dict{Symbol, Any}(Symbol(name) => Dict{String, Any}() for name in fieldnames(sh))
    rh = Dict{Symbol, Any}(Symbol(name) => Dict{String, Any}() for name in fieldnames(rh))
    ev = Dict{Symbol, Any}(Symbol(name) => Dict{String, Any}() for name in fieldnames(ev))
    ph = Dict{Symbol, Any}(Symbol(name) => Dict{String, Any}() for name in fieldnames(ph))
    d = Dict{Symbol, Any}(Symbol(name) => Dict{String, Any}() for name in fieldnames(d))
    
    bcs["T"] = [sh["T_s_in"][tt], ev["T_s_in"][tt], ph["T_s_in"][tt]]
    
    # Constants
    cp_const = 1443
    cp_slope = 0.172
    
    #---------------------------------------------------------
    bcs["cp_s"] = cp_const .+ cp_slope .* bcs["T"]  # [J/kg-K]
    #---------------------------------------------------------
    
    rho_salt = 2090 .- 0.636 .* T  # [kg/m3]
    cp_salt = 1443 .+ 0.172 .* T  # [J/kg-K]
    k_salt = 0.443 .+ 1.9e-4 .* T  # [W/m-K]
    mu_salt = (22.714 .- 0.120 .* T .+ 2.281e-4 .* T.^2 .- 1.474e-7 .* T.^3) ./ 1000  # [kg/m-s]
    Pr_salt = (mu_salt .* cp_salt) ./ k_salt
    
    # Discretize
    N_s1 = n["N_p"]
    N_s2 = N_s1 + n["N_sh"]
    N_s3 = N_s2 + n["N_rh"]
    N_s4 = N_s3 + n["N_ev"]
    N_s5 = N_s4 + n["N_ph"]
    
    # p
    p["rho_s"] = rho_salt[1:N_s1]
    p["cp_s"] = cp_salt[1:N_s1]
    p["k_s"] = k_salt[1:N_s1]
    p["mu_s"] = mu_salt[1:N_s1]
    p["Pr_s"] = Pr_salt[1:N_s1]
    
    # superheater
    sh["rho_s"] = rho_salt[N_s1+1:N_s2]
    sh["cp_s"] = cp_salt[N_s1+1:N_s2]
    sh["mu_s"] = mu_salt[N_s1+1:N_s2]
    sh["Pr_s"] = Pr_salt[N_s1+1:N_s2]
    sh["k_s"] = k_salt[N_s1+1:N_s2]
    
    # reheater
    rh["rho_s"] = rho_salt[N_s2+1:N_s3]
    rh["cp_s"] = cp_salt[N_s2+1:N_s3]
    rh["mu_s"] = mu_salt[N_s2+1:N_s3]
    rh["Pr_s"] = Pr_salt[N_s2+1:N_s3]
    rh["k_s"] = k_salt[N_s2+1:N_s3]
    
    # evaporator
    ev["rho_s"] = rho_salt[N_s3+1:N_s4]
    ev["cp_s"] = cp_salt[N_s3+1:N_s4]
    ev["mu_s"] = mu_salt[N_s3+1:N_s4]
    ev["Pr_s"] = Pr_salt[N_s3+1:N_s4]
    ev["k_s"] = k_salt[N_s3+1:N_s4]
    
    # preheater
    ph["rho_s"] = rho_salt[N_s4+1:N_s5]
    ph["cp_s"] = cp_salt[N_s4+1:N_s5]
    ph["mu_s"] = mu_salt[N_s4+1:N_s5]
    ph["Pr_s"] = Pr_salt[N_s4+1:N_s5]
    ph["k_s"] = k_salt[N_s4+1:N_s5]

    # Pipe 
    if any(p["Re_s"] .< 3000)
        p["Nu_s"] = 1.86 * (p["Re_s"] * p["Pr_s"] * p["D_i_p"] / p["dx"])^(1/3)
    elseif any((p["Re_s"] .>= 3000) .& (p["Re_s"] .<= 5E6))
        f = (100 * p["Re_s"])^(-1/4)   # blasius approx for hydraulically smooth ps
        p["Nu_s"] = (f / 8 * (p["Re_s"] - 1000) * p["Pr_s"]) / (1 + 12.7 * (f / 8)^0.5 * (p["Pr_s"]^(2/3) - 1))
    elseif any(p["Re_s"] > 5E6)
        f = (1.58 * log(p["Re_s"]) - 3.28)^(-2)
        p["Nu_s"] = (f / 2) * p["Re_s"] * p["Pr_s"] / (1.07 + 12.7 * (f / 2)^0.5 * (p["Pr_s"]^(2/3) - 1))
    end

    p["h_s"] = p["Nu_s"] * p["k_s"] / p["D_i_p"]  # Check Reynolds number range!!!!!! Re > 10,000

    ## Superheater
    function bell_delaware_method(sh::Dict)
        S_m     = sh["p_baffle"]*(sh["D_i_s"] - sh["D_otl"] + (sh["D_otl"] - sh["D_o_t"])*(sh["p_tube"] - sh["D_o_t"])/sh["p_tube"]);
        G_S     = sh["mdot_s(tt)"]/S_m;	# shell-side cross flow bulk velocity	
        sh["Re_s"]    = (sh["D_o_t"]*G_S)./sh["mu_s"];
        if  any(sh["Re_s < 10"])
            a_1     = 0.97;
            a_2     = -0.667;
            a_3     = 1.187;
            a_4     = 0.370;
        elseif any(sh["Re_s"] >= 10) .& any(sh["Re_s"] < 100)
            a_1     = 0.9;
            a_2     = -0.631;
            a_3     = 1.187;
            a_4     = 0.370;
        elseif any(sh["Re_s"] >= 100) .& any(sh["Re_s"] < 1000)
            a_1     = 0.408;
            a_2     = -0.460;
            a_3     = 1.187;
            a_4     = 0.370;
        elseif any(sh["Re_s"] >= 1000) .& any(sh["Re_s"] < 10000)
            a_1     = 0.107;
            a_2     = -0.266;
            a_3     = 1.187;
            a_4     = 0.370;
        elseif any(sh["Re_s"] >= 10000) 
            a_1     = 0.37;
            a_2     = -0.395;
            a_3     = 1.187;
            a_4     = 0.370; 
        end
        a           = a_3./(1+0.14.*sh["Re_s"].^a_4);
        j           = a_1*((1.33)/(sh["p_tube"]/(sh["D_o_t"]))).^a.*sh["Re_s"].^a_2;
        mu_wall 	= sh["mu_s"];
        h_ideal 	= j.*sh["cp_s"].*G_S.*sh["Pr_s"].^(-2/3)...
                        .*(sh["mu_s"])./mu_wall;
        P_T_eff     = sh["p_tube"];
        D_ctl       = sh["D_otl"] - sh["D_o_t"];
        theta_ctl	= 2*acos(sh["D_i_s"]*(1 - 2*sh["cut"])/D_ctl);
        F_W         = (theta_ctl - sin(theta_ctl))/(2*pi);
        F_C         = 1 - 2*F_W;
        J_C         = 0.55 + 0.72*F_C;
        theta_ds    = 2*acos(1 - 2*sh["cut"]);
        S_sb        = sh["D_i_s"]*(sh["D_i_s"] - sh["D_baffle"])*(pi -0.5*theta_ds);
        S_tb        = (pi/4)*(sh["D_bch"]^2 - sh["D_o_t"]^2)*sh["N_tubes"]*(1-F_W);
        ratio_S     = S_sb / (S_sb + S_tb);
        ratio_L     = (S_sb + S_tb)/S_m;
        J_L         = 0.44*(1-ratio_S) + (1-0.44*(1-ratio_S))*exp(-2.2*ratio_L);
        L_TB        = sh["D_bch"] - sh["D_o_t"];
        C_j         = 1.25;
        S_b         = sh["p_baffle"]*(sh["D_i_s"] - sh["D_otl"] - sh["D_o_t"]/2); # bundle by-pass flow
        ratio_ss    = 0;
        n1          = 0.6; 	# Re >= 100
        J_B         = exp(-C_j*(S_b/S_m)*(1 - (2*ratio_ss)^(1/3)));
        N_tcc       = (sh["D_i_s"]/sh["p_tube"])*(1 - 2*sh["cut"]);
        J_R         = 1;	# adverse temp gradiant
        N_tcw       = (0.8/sh["p_tube"])*(sh["D_i_s"]*sh["cut"] - (sh["D_i_s"] - (sh["D_otl"] - sh["D_o_t"]))/2);
        J_S         = ((sh["N_baffles"] - 1)+(sh["p_baffle"]/sh["p_baffle"])^(1-n1) ...
                            + (sh["p_baffle"]/sh["p_baffle"])^(1-n1)) ...
                        /((sh["N_baffles"]-1) + (sh["p_baffle"]/sh["p_baffle"]) ...
                            + (sh["p_baffle"]/sh["p_baffle"])); #"assumed baffle pitch in/out = paffle pitch"
        sh["h_s"]      = h_ideal.*J_C.*J_R.*J_L.*J_B.*J_S; 
        return sh 
    end 
    sh = bell_delaware_method(sh)
    rh = bell_delaware_method(rh)
    ev = bell_delaware_method(ev)
    ph = bell_delaware_method(ph)
    # Calculate pipe parameters for each section
    # Note: Re calculations and other conditional statements remain unchanged
    
    return bcs, p, sh, rh, ev, ph
end
