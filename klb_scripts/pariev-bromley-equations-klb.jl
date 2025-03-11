# Following equations are taken from Pariev & Bromley 1998 prescription of velocity. 

# Eddington Luminosity function
function eddington_luminosity(M)
    Ledd_constant = 1.2e46  # Base Eddington luminosity in ergs/s for 10^8 solar masses
    return Ledd_constant * (M / 1e8)  # M should be in units of solar masses (M☉)
end

# Equation 4
function A(r_star, a_star)
    return 1 + ((a_star^2) / r_star^2) + (2 * (a_star^2) / r_star^3) 
end

# Equation 5
function B(r_star, a_star)
    return 1 + (a_star / r_star^(3/2))
end

# Equation 6
function C(r_star, a_star)
    return 1 - (3 / r_star) + ((2 * a_star) / r_star^(3/2))
end

# Equation 7
function D(r_star, a_star)
    return 1 - (2 / r_star) + ((a_star^2) / r_star^2)
end

# Equation 8
function E(r_star, a_star)
    return 1 + (4 * (a_star^2) / r_star^2) - (4 * (a_star^2) / r_star^3) + (3 * (a_star^4) / r_star^4)
end

# Equation 9
function F(r_star, a_star)
    return 1 - (2 * a_star / r_star^(3/2)) + (a_star^2 / r_star^2)
end

# Equation 10
function G(r_star, a_star)
    return 1 - (2 / r_star) + (a_star / r_star^(3/2))
end

# Equation 11
function J(r_star, a_star)
    numerator = 1 + (a_star / r_star^(3/2))
    denominator = (1 - (3 / r_star) + (2 * a_star / r_star^(3/2)))^(1/2)
    return numerator / denominator
end

# Equation 12
function L_func(r_star, a_star, r_ms) # ('_func' to not mistake with L, the luminosity)
    term1 = F(r_star, a_star) / sqrt(C(r_star, a_star))
    inner_term = 1 - (2 * a) / (3 * sqrt(r_ms))
    if inner_term < 0
        error("The argument inside the root becomes negative.")
    end
    term2 = (2 * sqrt(3) / sqrt(r_star)) * inner_term
    return term1 - term2
end

# Equation 35 from Page & Thorne (1975)
# M is geometrised units, (i.e.= 1, for the sake of calculating the ISCO)
function Q(r, a_star, M)
    # Validate inputs
    if r <= 0 || M <= 0
        error("Radius (r) and mass (M) must be positive.")
    end

    # Clamp 'a' to valid range for acos
    a_star = clamp(a_star, -1.0, 1.0)

    R_isco = Gradus.isco(KerrMetric(M, a_star))
    x = sqrt(r / M)
    x0 = sqrt(R_isco / M)
    x1 = 2 * cos((1 / 3) * (acos(a_star)) - (π / 3))
    x2 = 2 * cos((1 / 3) * (acos(a_star)) + (π / 3))
    x3 = -2 * cos((1 / 3) * (acos(a_star)))

    # Ensure log arguments are positive
    if x <= x1 || x <= x2 || x <= x3 || x0 <= x1 || x0 <= x2 || x0 <= x3
        error("Logarithm argument becomes non-positive. Check inputs.")
    end

    numerator = 1 + a_star * x^(-3)
    denominator = sqrt(1 - 3 * x^(-2) + 2 * a_star * x^(-3))
    prefactor = numerator / denominator * (1 / x)
    
    term1 = (x - x0) - (3/2) * a_star * log(x / x0)
    term2 = -(3 * (x1 - a_star)^2) / (x1 * (x1 - x2) * (x1 - x3)) * log((x - x1) / (x0 - x1))
    term3 = -(3 * (x2 - a_star)^2) / (x2 * (x2 - x1) * (x2 - x3)) * log((x - x2) / (x0 - x2))
    term4 = -(3 * (x3 - a_star)^2) / (x3 * (x3 - x1) * (x3 - x2)) * log((x - x3) / (x0 - x3))


    Q_val = prefactor * (term1 + term2 + term3 + term4)

    return Q_val
end


# Equation 15: Speed of sound to speed of light ratio (v_turb)
function sound_speed_ratio(r, a, epsilon, L, L_edd, r_ms, M)
    if r < r_ms
        return 0
    end

    r_star = r/M 
    a_star= a/M 

    factor = 1.18 * epsilon^(-1) * (L / L_edd) * r^(-3/2)
    A_val = A(r_star, a_star)
    B_val = B(r_star, a_star)
    D_val = D(r_star, a_star)
    E_val = E(r_star, a_star)
    Q_val = Q(r, a_star, M)

    return (factor * A_val * B_val^(-2) * D_val^(-1/2) * E_val^(-1/2) * Q_val) 
end
