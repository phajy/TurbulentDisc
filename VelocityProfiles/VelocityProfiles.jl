using Gradus
using Plots

function A(r, a, M)

    r_star = r/M
    a_star = a/M

    return 1 + (a_star^2)/(r_star^2) + (2*(a_star^2))/(r_star^3)

end

function B(r, a, M)

    r_star = r/M
    a_star = a/M
    
    return 1 + (a_star)/(r_star^(3/2))

end

function C(r, a, M)

    r_star = r/M
    a_star = a/M
    
    return 1 - (3)/(r_star) + (2*a_star)/(r_star^(3/2))

end

function D(r, a, M)

    r_star = r/M
    a_star = a/M
    
    return 1 - (2)/(r_star) + (a_star^2)/(r_star^2)

end

function E(r, a, M)

    r_star = r/M
    a_star = a/M
    
    return 1 + (4*(a_star^2))/(r_star^2) - (4*(a_star^2))/(r_star^3) + (3*(a_star^4))/(r_star^4)

end

function F(r, a, M)

    r_star = r/M
    a_star = a/M
    
    return 1 - (2*(a_star))/(r_star^(3/2)) + (a_star^2)/(r_star^2)

end

function G(r, a, M)

    r_star = r/M
    a_star = a/M
    
    return 1 - (2)/(r_star) + (a_star^2)/(r_star^(3/2))

end

function I(r, a, M)

    r_star = r/M
    a_star = a/M
    
    return (1 + (a_star)/(r_star^(3/2)))/((1 - 3/r_star + (2*a_star)/(r_star^(3/2)))^(1/2))

end

function L(r, a, M)

    r_star = r/M
    a_star = a/M
    
    return F(r, a, M)/(C(r, a, M)^(1/2)) - (2*sqrt(3))/(sqrt(r_star))*(1 - (2*a_star)/(3*sqrt(r_star)))
    
end

function Ledd(M)

    M_sol = 1

    return  1.2*(10^46) * ((M)/(M_sol*(10^8)))

end

# Q(r, a, M) as defined in Page & Thorne (1974), their equation 35
function Q(r, a, M, isco)
    return L(r, a, M)
    # note above line is wrong!
end

function SoundSpeed(r, a, M, ϵ, isco)

    r_star = r/M
    a_star = a/M

    L_Ledd = 1
    LComp = 1

    return 1.18 * (ϵ^(-1)) * r_star^(-3/2) * A(r, a, M) * B(r, a, M)^(-2) * D(r, a, M)^(-1/2) * E(r, a, M)^(-1/2) * Q(r, a, M, isco)

end

function RadialSpeed(r, a, alpha, M, ϵ, isco)

    r_star = r/M
    a_star = a/M

    L_Ledd = 1
    LComp = 1

    return 1.18 * alpha * (ϵ^(-2)) * r_star^(-5/2) * A(r, a, M)^2 * B(r, a, M)^(-3) * C(r, a, M)^(-1/2) * D(r, a, M)^(-1/2) * E(r, a, M)^(-1) * Q(r, a, M, isco)

end

a = 0.9
M = 1.0
alpha = 0.3

# calculate the ISCO and efficiency for the given spin parameter
# note that ISCO and r_ms are used interchangeably in the literature
m = KerrMetric(M = M, a = a)
isco = Gradus.isco(m)
ϵ = 1.0 - Gradus.CircularOrbits.energy(m, isco)

rPos = collect(range(isco, 100, 100))

cs = [SoundSpeed(r, a, M, ϵ, isco) for r in rPos]
vr = [RadialSpeed(r, a, alpha, M, ϵ, isco) for r in rPos]

plot(rPos, cs)
# plot(rPos, vr)