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

# Equation 35 from Page & Thorne (1975)
# M is geometrised units, (i.e.= 1, for the sake of calculating the ISCO)
function Q(r, a, M)
    # Validate inputs
    if r <= 0 || M <= 0
        error("Radius (r) and mass (M) must be positive.")
    end

    # Clamp 'a' to valid range for acos
    a = clamp(a, -1.0, 1.0)

    R_isco = Gradus.isco(KerrMetric(M, a))
    x = sqrt(r / M)
    x0 = sqrt(R_isco / M)
    x1 = 2 * cos((1 / 3) * (acos(a)) - (π / 3))
    x2 = 2 * cos((1 / 3) * (acos(a)) + (π / 3))
    x3 = -2 * cos((1 / 3) * (acos(a)))

    # Ensure log arguments are positive
    if x <= x1 || x <= x2 || x <= x3 || x0 <= x1 || x0 <= x2 || x0 <= x3
        error("Logarithm argument becomes non-positive. Check inputs.")
    end

    # Ensure denominators are non-zero
    denom1 = x1 * (x1 - x2) * (x1 - x3)
    denom2 = x2 * (x2 - x1) * (x2 - x3)
    denom3 = x3 * (x3 - x1) * (x3 - x2)
    if denom1 == 0 || denom2 == 0 || denom3 == 0
        error("Denominator in Q(r, a, M) becomes zero. Check inputs.")
    end

    Q_val =
        (3 / (2 * M)) *
        (1 / (x^2 * (x^3 - (3 * x) + (2 * a)))) *
        (
            x - x0 - ((3 / 2) * a * log(x / x0)) -
            (3 * (x1 - a)^2) / denom1 * log((x - x1) / (x0 - x1)) -
            (3 * (x2 - a)^2) / denom2 * log((x - x2) / (x0 - x2)) -
            (3 * (x3 - a)^2) / denom3 * log((x - x3) / (x0 - x3))
        )
    return Q_val
end



# Equation 15: Speed of sound to speed of light ratio (v_turb)
function sound_speed_ratio(r, a, epsilon, L, L_edd, r_ms, M)
    if r < r_ms
        return 0
    end
    factor = 1.18 * epsilon^(-1) * (L / L_edd) * r^(-3/2)
    A_val = A(r, a)
    B_val = B(r, a)
    D_val = D(r, a)
    E_val = E(r, a)
    Q_val = Q(r, a, M)

    return factor * A_val * B_val^(-2) * D_val^(-1/2) * E_val^(-1/2) * Q_val
end
