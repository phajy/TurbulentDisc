"""
This script defines the sound speed function (Equation 15) and radial inflow velocity function (Equation 16),
as described by Pariev & Bromley (1998).
"""
using Gradus

function A(r_star, a_star, M)

    return 1 + (a_star^2)/(r_star^2) + (2*(a_star^2))/(r_star^3)

end

function B(r_star, a_star, M)
    
    return 1 + (a_star)/(r_star^(3/2))

end

function C(r_star, a_star, M)
    
    return 1 - (3)/(r_star) + (2*a_star)/(r_star^(3/2))

end

function D(r_star, a_star, M)
    
    return 1 - (2)/(r_star) + (a_star^2)/(r_star^2)

end

function E(r_star, a_star, M)
    
    return 1 + (4*(a_star^2))/(r_star^2) - (4*(a_star^2))/(r_star^3) + (3*(a_star^4))/(r_star^4)

end

function F(r_star, a_star, M)
    
    return 1 - (2*(a_star))/(r_star^(3/2)) + (a_star^2)/(r_star^2)

end

function G(r_star, a_star, M)
    
    return 1 - (2)/(r_star) + (a_star)/(r_star^(3/2))

end

function I(r_star, a_star, M)
    
    return (1 + (a_star)/(r_star^(3/2))) / ( ( 1 - 3/r_star + (2*a_star)/(r_star^(3/2)) ) ^(1/2) )

end

function L(r_star, a_star, M)
    
    return F(r_star, a_star, M)/(C(r_star, a_star, M)^(1/2)) - ( (2*sqrt(3)) / (sqrt(r_star)) ) * (1 - (2*a_star)/(3*sqrt(r_star)))
    
end

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

function efficiency(m)

    return 1 - Gradus.CircularOrbits.energy(m, Gradus.isco(m))

end

function SoundSpeed(m, r, a, M, ratio)

    if r < Gradus.isco(m)
        return 0
    end

    r_star = r/M
    a_star = a/M

    return 1.18 *
    (efficiency(m)^(-1)) *
    (ratio) *
    r_star^(-3/2) *
    A(r_star, a_star, M) *
    B(r_star, a_star, M)^(-2) *
    D(r_star, a_star, M)^(-1/2) *
    E(r_star, a_star, M)^(-1/2) *
    Q(r, a_star, M)

end

function RadialSpeed(m, r, a, alpha, M, ratio)

    if r < Gradus.isco(m)
        return 0
    end

    r_star = r/M
    a_star = a/M

    return 1.13 *
    alpha *
    (efficiency(m)^(-2)) *
    (ratio)^2 *
    r_star^(-5/2) *
    A(r_star, a_star, M)^2 *
    B(r_star, a_star, M)^(-3) *
    C(r_star, a_star, M)^(-1/2) *
    D(r_star, a_star, M)^(-1/2) *
    E(r_star, a_star, M)^(-1) *
    Q(r, a, M)

end