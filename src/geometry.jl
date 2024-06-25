using MAT
using LinearAlgebra

function geometry(d)
    G = matread("G.mat")["G"]
    
    p = Dict{Symbol, Any}()
    sh = Dict{Symbol, Any}()
    rh = Dict{Symbol, Any}()
    ev = Dict{Symbol, Any}()
    ph = Dict{Symbol, Any}()

    # p geometry
    p[:L] = 30.0  # [m] length from measurement location to superheater inlet
    p[:th_p] = 0.034  # [m] p thickness
    p[:D_i_p] = 0.234  # [m] inner diameter of p
    p[:D_o_p] = p[:D_i_p] + 2 * p[:th_p]  # [m] outer diameter of p
    p[:r_i_p] = p[:D_i_p] / 2  # [m] inner radius of p
    p[:r_c_p] = p[:r_i_p] + (p[:th_p] / 2)
    p[:r_o_p] = p[:D_o_p] / 2  # [m] outer radius of p
    p[:A_c_mw] = π * (p[:r_o_p]^2 - p[:r_i_p]^2)  # [m^2] cross-sectional area of metal
    p[:A_c_sa] = π * p[:r_i_p]^2  # [m^2] cross-sectional area of salt
    p[:th_insul] = 0.060  # insulation thickness
    p[:D_o_insul] = p[:D_o_p] + 2 * p[:th_insul]
    p[:r_o_insul] = p[:D_o_insul] / 2

    # Superheater geometry
    sh[:N_tubes] = G[18, 4]             # [#] of tubes
    sh[:N_baffles] = G[34, 4]           # [#] of baffles
    sh[:L] = 2 * G[14, 4] + d[:elbow]   # [m] average tube length
    sh[:D_o_t] = G[4, 4]                # [m] outer diameter of tube
    sh[:D_i_t] = 2 * G[1, 4]            # [m] inner diameter of tube
    sh[:D_i_s] = 2 * G[19, 4]           # [m] inner shell diameter
    sh[:D_otl] = G[17, 4]               # [m] tube bundle diameter
    sh[:D_baffle] = G[33, 4]  # [m] baffle diameter
    sh[:layout] = G[16, 4]  # layout of tube pitch
    sh[:p_tube] = G[15, 4]  # [m] tube pitch
    sh[:p_baffle] = G[27, 4]  # [m] baffle pitch
    sh[:cut] = G[30, 4]  # [%] baffle cut
    sh[:D_bch] = G[32, 4]  # [m] baffle clearance hole diameter
    sh[:P_s] = 2 * π * G[2, 4]  # [m] perimeter of salt
    sh[:P_w] = 2 * π * G[1, 4]  # [m] perimeter of water
    sh[:A_c_s] = π * G[19, 4]^2 - π * G[2, 4]^2 * sh[:N_tubes]  # [m^2] cross-sectional area of salt
    sh[:A_c_m] = π * (G[2, 4]^2 - G[1, 4]^2) * sh[:N_tubes]  # [m^2] cross-sectional area of metal
    sh[:A_c_w] = π * G[1, 4]^2 * sh[:N_tubes]  # [m^2] cross-sectional area of water
    sh[:V_h] = π * G[19, 4]^2 * 0.19  # [m^2] Volume of header/tubesheet
    sh[:A_c_shell] = π * (G[20, 4]^2 - G[19, 4]^2)  # [m^2] cross-sectional area of metal shell
    sh[:A_c_shell_i] = G[25, 4]
    sh[:r_tube_o] = G[2, 4]
    sh[:r_tube_i] = G[1, 4]
    sh[:r_shell_o] = G[20, 4]
    sh[:r_shell_i] = G[19, 4]
    sh[:th_tube] = sh[:r_tube_o] - sh[:r_tube_i]
    sh[:r_tube_c] = sh[:r_tube_o] - sh[:th_tube] / 2
    sh[:A_c_tube] = π * sh[:r_tube_i]^2
    sh[:th_header] = 0.016
    sh[:th_tbs] = 0.19
    sh[:h] = (sh[:p_tube] - sh[:D_o_t]) / 2
    sh[:R] = sh[:p_tube] / 2
    sh[:nu_lig] = sh[:h] / sh[:R]  # ligament efficiency
    sh[:Npass] = 1  # number of passes

    # Reheater geometry
    rh[:N_tubes] = G[18, 3]  # [#] of tubes
    rh[:N_baffles] = G[34, 3]  # [#] of baffles
    rh[:L] = 2 * G[14, 3] + 1  # [m] average tube length
    rh[:D_o_t] = G[4, 3]  # [m] outer diameter of tube
    rh[:D_i_t] = 2 * G[1, 3]  # [m] outer diameter of tube
    rh[:D_i_s] = 2 * G[19, 3]  # [m] inner shell diameter
    rh[:D_otl] = G[17, 3]  # [m] tube bundle diameter
    rh[:D_baffle] = G[33, 3]  # [m] baffle diameter
    rh[:layout] = G[16, 3]  # layout of tube pitch
    rh[:p_tube] = G[15, 3]  # [m] tube pitch
    rh[:p_baffle] = G[27, 3]  # [m] baffle pitch
    rh[:cut] = G[30, 3]  # [%] baffle cut
    rh[:D_bch] = G[32, 3]  # [m] baffle clearance hole diameter
    rh[:P_s] = 2 * π * G[2, 3]  # [m] perimeter of salt
    rh[:P_w] = 2 * π * G[1, 3]  # [m] perimeter of water
    rh[:A_c_s] = π * G[19, 3]^2 - π * G[2, 3]^2 * rh[:N_tubes]  # [m^2] cross-sectional area of salt
    rh[:A_c_m] = π * (G[2, 3]^2 - G[1, 3]^2) * rh[:N_tubes]  # [m^2] cross-sectional area of metal
    rh[:A_c_w] = π * G[1, 3]^2 * rh[:N_tubes]  # [m^2] cross-sectional area of water
    rh[:V_h] = π * G[19, 3]^2 * 0.19  # [m^2] Volume of header/tubesheet
    rh[:A_c_shell] = π * (G[20, 3]^2 - G[19, 3]^2)  # [m^2] cross-sectional area of metal shell
    rh[:A_c_shell_i] = G[25, 3]
    rh[:r_tube_o] = G[2, 3]
    rh[:r_tube_i] = G[1, 3]
    rh[:r_shell_o] = G[20, 3]
    rh[:r_shell_i] = G[19, 3]
    rh[:th_tube] = rh[:r_tube_o] - rh[:r_tube_i]
    rh[:r_tube_c] = rh[:r_tube_o] - rh[:th_tube] / 2
    rh[:A_c_tube] = π * rh[:r_tube_i]^2
    rh[:h] = (rh[:p_tube] - rh[:D_o_t]) / 2
    rh[:R] = rh[:p_tube] / 2
    rh[:th_tbs] = 0.19
    rh[:Lig_eff] = rh[:h] / rh[:R]
    rh[:Npass] = 1  # number of passes

    # Evaporator geometry
    ev[:N_tubes] = G[18, 2]  # [#] of tubes
    ev[:N_baffles] = G[34, 2]  # [#] of baffles
    ev[:L] = 2 * G[14, 2]  # [m] average tube length
    ev[:D_o_t] = G[4, 2]  # [m] outer diameter of tube
    ev[:D_i_t] = 2 * G[1, 2]  # [m] inner diameter of tube
    ev[:D_i_s] = 2 * G[19, 2]  # [m] inner shell diameter
    ev[:D_otl] = G[17, 2]  # [m] tube bundle diameter
    ev[:D_baffle] = G[33, 2]  # [m] baffle diameter
    ev[:layout] = G[16, 2]  # layout of tube pitch
    ev[:p_tube] = G[15, 2]  # [m] tube pitch
    ev[:p_baffle] = G[27, 2]  # [m] baffle pitch
    ev[:cut] = G[30, 2]  # [%] baffle cut
    ev[:D_bch] = G[32, 2]  # [m] baffle clearance hole diameter
    ev[:P_s] = 2 * π * G[2, 2]  # [m] perimeter of salt
    ev[:P_w] = 2 * π * G[1, 2]  # [m] perimeter of water
    ev[:A_c_s] = π * G[19, 2]^2 - π * G[2, 2]^2 * ev[:N_tubes]  # [m^2] cross-sectional area of salt
    ev[:A_c_m] = π * (G[2, 2]^2 - G[1, 2]^2) * ev[:N_tubes]  # [m^2] cross-sectional area of metal
    ev[:A_c_w] = π * G[1, 2]^2 * ev[:N_tubes]  # [m^2] cross-sectional area of water
    ev[:V_h] = π * G[19, 2]^2 * 0.19  # [m^2] Volume of header/tubesheet
    ev[:A_c_shell] = π * (G[20, 2]^2 - G[19, 2]^2)  # [m^2] cross-sectional area of metal shell
    ev[:A_c_shell_i] = G[25, 2]
    ev[:r_tube_o] = G[2, 2]
    ev[:r_tube_i] = G[1, 2]
    ev[:r_shell_o] = G[20, 2]
    ev[:r_shell_i] = G[19, 2]
    ev[:th_tube] = ev[:r_tube_o] - ev[:r_tube_i]
    ev[:r_tube_c] = ev[:r_tube_o] - ev[:th_tube] / 2
    ev[:A_c_tube] = π * ev[:r_tube_i]^2
    ev[:h] = (ev[:p_tube] - ev[:D_o_t]) / 2
    ev[:R] = ev[:p_tube] / 2
    ev[:th_tbs] = 0.19
    ev[:Lig_eff] = ev[:h] / ev[:R]
    ev[:lane] = 0.048 * 2  # diametral lane distance (on tubesheet)
    ev[:K_D] = ev[:r_tube_i] / ev[:lane]
    ev[:Npass] = 1  # number of passes

    # Preheater geometry
    ph[:N_tubes] = G[18, 1]  # [#] of tubes
    ph[:N_baffles] = G[34, 1]  # [#] of baffles
    ph[:L] = 2 * G[14, 1]  # [m] average tube length
    ph[:D_o_t] = G[4, 1]  # [m] outer diameter of tube
    ph[:D_i_t] = 2 * G[1, 1]  # [m] inner diameter of tube
    ph[:D_i_s] = 2 * G[19, 1]  # [m] inner shell diameter
    ph[:D_otl] = G[17, 1]  # [m] tube bundle diameter
    ph[:D_baffle] = G[33, 1]  # [m] baffle diameter
    ph[:layout] = G[16, 1]  # layout of tube pitch
    ph[:p_tube] = G[15, 1]  # [m] tube pitch
    ph[:p_baffle] = G[27, 1]  # [m] baffle pitch
    ph[:cut] = G[30, 1]  # [%] baffle cut
    ph[:D_bch] = G[32, 1]  # [m] baffle clearance hole diameter
    ph[:P_s] = 2 * π * G[2, 1]  # [m] perimeter of salt
    ph[:P_w] = 2 * π * G[1, 1]  # [m] perimeter of water
    ph[:A_c_s] = π * G[19, 1]^2 - π * G[2, 1]^2 * ph[:N_tubes]  # [m^2] cross-sectional area of salt
    ph[:A_c_m] = π * (G[2, 1]^2 - G[1, 1]^2) * ph[:N_tubes]  # [m^2] cross-sectional area of metal
    ph[:A_c_w] = π * G[1, 1]^2 * ph[:N_tubes]  # [m^2] cross-sectional area of water
    ph[:V_h] = π * G[19, 1]^2 * 0.19  # [m^2] Volume of header/tubesheet
    ph[:A_c_shell] = π * (G[20, 1]^2 - G[19, 1]^2)  # [m^2] cross-sectional area of metal shell
    ph[:A_c_shell_i] = G[25, 1]
    ph[:r_tube_o] = G[2, 1]
    ph[:r_tube_i] = G[1, 1]
    ph[:r_shell_o] = G[20, 1]
    ph[:r_shell_i] = G[19, 1]
    ph[:th_tube] = ph[:r_tube_o] - ph[:r_tube_i]
    ph[:r_tube_c] = ph[:r_tube_o] - ph[:th_tube] / 2
    ph[:A_c_tube] = π * ph[:r_tube_i]^2
    ph[:h] = (ph[:p_tube] - ph[:D_o_t]) / 2
    ph[:R] = ph[:p_tube] / 2
    ph[:th_tbs] = 0.19
    ph[:Lig_eff] = ph[:h] / ph[:R]
    ph[:lane] = 0.096  # diametral lane distance (on tubesheet)
    ph[:K_D] = ph[:r_tube_i] / ph[:lane]
    ph[:Npass] = 1  # number of passes

    return p, sh, rh, ev, ph, d
end

# Example of calling the function with a sample input
d = Dict(:elbow => 0.5)  # replace with actual data
p, sh, rh, ev, ph, d = geometry(d)