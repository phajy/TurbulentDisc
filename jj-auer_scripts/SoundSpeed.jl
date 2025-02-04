using Plots
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
    
    return (1 + (a_star)/(r_star^(3/2)))/((1 - 3/r_star + (2*a_star)/(r_star^(3/2)))^(1/2))

end

function L(r_star, a_star, M)
    
    return F(r_star, a_star, M)/(C(r_star, a_star, M)^(1/2)) - ((2*sqrt(3))/(sqrt(r_star)))*(1 - (2*a_star)/(3*sqrt(r_star)))
    
end

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
        if x <= x1
            error("x <= x1")
        elseif x <= x2
            error("x <= x2")
        elseif x <= x3
            error("x <= x3")
        elseif x0 <= x1
            error("x0 <= x1")
        elseif x0 <= x2
            error("x0 <= x2")
        elseif x0 <= x3
            error("x0 <= x3")
        end

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

function efficiency(m)

    return 1 - Gradus.CircularOrbits.energy(m, Gradus.isco(m))

end

function LumEdd(M)

    return  1.2e46 * (M)

end



function SoundSpeed(m, r, a, M, lum)

    if r < Gradus.isco(m)
        return 0
    end

    r_star = r/M
    a_star = a/M

    return 1.18 *
    (efficiency(m)^(-1)) *
    (lum/LumEdd(M)) *
    r_star^(-3/2) *
    A(r_star, a_star, M) *
    B(r_star, a_star, M)^(-2) *
    D(r_star, a_star, M)^(-1/2) *
    E(r_star, a_star, M)^(-1/2) *
    Q(r, a, M)

end

function RadialSpeed(m, r, a, alpha, M, lum)

    if r < Gradus.isco(m)
        return 0
    end

    r_star = r/M
    a_star = a/M

    return 1.13 *
    alpha *
    (efficiency(m)^(-2)) *
    (Lum/LumEdd(M))^2 *
    r_star^(-5/2) *
    A(r_star, a_star, M)^2 *
    B(r_star, a_star, M)^(-3) *
    C(r_star, a_star, M)^(-1/2) *
    D(r_star, a_star, M)^(-1/2) *
    E(r_star, a_star, M)^(-1) *
    Q(r, a, M)

end

M = 1.0
alpha = 0.1
lum = 1e46

plt = plot(
    xlabel = "Radius (r/M)",
    ylabel = "Sound Speed (v/c)",
    title = "Sound Speed vs. Radius",
    legend = :topright,
)

rPos = collect(range(1.0, stop=25.0, length=500))

for a in [0, 0.5, 0.90, 0.99, 0.998]
    m = KerrMetric(M = M, a = a)
    plot!(rPos, [SoundSpeed(m, r, a, M, lum) for r in rPos], label="a/M = $a")
end

display(plt)
