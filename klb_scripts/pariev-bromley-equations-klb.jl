# Following equations are taken from Pariev & Bromley 1998 prescription of velocity. 

# Eddington Luminosity function
function eddington_luminosity(M)
    Ledd_constant = 1.2e46  # Base Eddington luminosity in ergs/s for 10^8 solar masses
    return Ledd_constant * (M / 1e8)  # M should be in units of solar masses (M☉)
end

# Equation 4
function A(r, a)
    return 1 + ((a^2) / r^2) + (2 * (a^2) / r^3) 
end

# Equation 5
function B(r, a)
    return 1 + (a / r^(3/2))
end

# Equation 6
function C(r, a)
    return 1 - (3 / r) + ((2 * a) / r^(3/2))
end

# Equation 7
function D(r, a)
    return 1 - (2 / r) + ((a^2) / r^2)
end

# Equation 8
function E(r, a)
    return 1 + (4 * (a^2) / r^2) - (4 * (a^2) / r^3) + (3 * (a^4) / r^4)
end

# Equation 9
function F(r, a)
    return 1 - (2 * a / r^(3/2)) + (a^2 / r^2)
end

# Equation 10
function G(r, a)
    return 1 - (2 / r) + (a / r^(3/2))
end

# Equation 11
function J(r, a)
    numerator = 1 + (a / r^(3/2))
    denominator = (1 - (3 / r) + (2 * a / r^(3/2)))^(1/2)
    return numerator / denominator
end

# Equation 12
function L_func(r, a, r_ms) # ('_func' to not mistake with L, the luminosity)
    term1 = F(r, a) / sqrt(C(r, a))
    inner_term = 1 - (2 * a) / (3 * sqrt(r_ms))
    if inner_term < 0
        error("The argument inside the root becomes negative.")
    end
    term2 = (2 * sqrt(3) / sqrt(r)) * inner_term
    return term1 - term2
end

# Equation 15: Speed of sound to speed of light ratio (v_turb)
function sound_speed_ratio(r, a, epsilon, L, L_edd, r_ms)
    if r < r_ms
        return 0
    end
    factor = 1.18 * epsilon^(-1) * (L / L_edd) * r^(-3/2)
    A_val = A(r, a)
    B_val = B(r, a)
    D_val = D(r, a)
    E_val = E(r, a)
    L_val = L_func(r, a, r_ms)
    return factor * A_val * B_val^(-2) * D_val^(-1/2) * E_val^(-1/2) * L_val
end
