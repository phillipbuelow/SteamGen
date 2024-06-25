function insulation(T)
    k_insul = 0.0382178 + 0.0000431732 * T + 2.58496e-7 * T^2
    rho_insul = 128  # [kg/m3]
    return k_insul, rho_insul
end
