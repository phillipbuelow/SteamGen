@kwdef struct sh
    # loading example industry scale STHXs
    G = JLD2.load("data/geometry/geo_example.jld2")["G"]

    N_tubes::Int64    = G[18, 4]             # [#] of tubes
    N_baffles::Int64  = G[34, 4]           # [#] of baffles
    L::Float64          = 2 * G[14, 4] + d[:elbow]   # [m] average tube length
    D_o_t::Float64      = G[4, 4]                # [m] outer diameter of tube
    D_i_t::Float64      = 2 * G[1, 4]            # [m] inner diameter of tube
    D_i_s::Float64      = 2 * G[19, 4]           # [m] inner shell diameter
    D_otl::Float64      = G[17, 4]               # [m] tube bundle diameter
    D_baffle::Float64   = G[33, 4]  # [m] baffle diameter
    layout::Float64     = G[16, 4]  # layout of tube pitch
    p_tube::Float64     = G[15, 4]  # [m] tube pitch
    p_baffle::Float64   = G[27, 4]  # [m] baffle pitch
    cut::Float64        = G[30, 4]  # [%] baffle cut
    D_bch::Float64      = G[32, 4]  # [m] baffle clearance hole diameter
    P_s::Float64        = 2 * π * G[2, 4]  # [m] perimeter of salt
    P_w::Float64        = 2 * π * G[1, 4]  # [m] perimeter of water
    A_c_s::Float64      = π * G[19, 4]^2 - π * G[2, 4]^2 * N_tubes  # [m^2] cross-sectional area of salt
    A_c_m::Float64      = π * (G[2, 4]^2 - G[1, 4]^2) * N_tubes  # [m^2] cross-sectional area of metal
    A_c_w::Float64      = π * G[1, 4]^2 * N_tubes  # [m^2] cross-sectional area of water
    V_h::Float64        = π * G[19, 4]^2 * 0.19  # [m^2] Volume of header/tubesheet
    A_c_shell::Float64  = π * (G[20, 4]^2 - G[19, 4]^2)  # [m^2] cross-sectional area of metal shell
    A_c_shell_i::Float64 = G[25, 4]
    r_tube_o::Float64   = G[2, 4]
    r_tube_i::Float64   = G[1, 4]
    r_shell_o::Float64  = G[20, 4]
    r_shell_i::Float64  = G[19, 4]
    th_tube::Float64    = r_tube_o - r_tube_i
    r_tube_c::Float64   = r_tube_o - th_tube / 2
    A_c_tube::Float64   = π * r_tube_i^2
    th_header::Float64  = 0.016
    th_tbs::Float64     = 0.19
    h::Float64          = (p_tube - D_o_t) / 2
    R::Float64          = p_tube / 2
    nu_lig::Float64     = h / R                     # ligament efficiency
    Npass::Int64        = 1                         # number of passes
end