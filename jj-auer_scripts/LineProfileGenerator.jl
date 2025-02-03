# Script to plot both the zero-turbulence thin disc and perlin/fBm noise turbulence line profiles for different emissivities 
# and inclination angles (emissivities and inclinations angle ranges as of Pariev & Bromley 1998)

# Import libraries
using Plots, Gradus
include("TurbulenceMaps.jl")

# Functions used for turbulence
# +------------------------------------------------------------------------------+

function turbulent_redshift(metric, x_obs, vel_func, a, M, lum, correlation_length, mach)
    # metric matrix at the observer's position
    g_obs = Gradus.metric(metric, x_obs)
    # fixed stationary observer velocity
    v_obs = SVector{4,eltype(x)}(1, 0, 0, 0)

    # internal closure
    function _internal_turbulent_redshift(m::AbstractMetric, gp, t)
        v_disc = vel_func(m, gp.x[2], gp.x[4], a, M, lum, correlation_length, mach)

        g = Gradus.metric(m, gp.x)
        Gradus.RedshiftFunctions._redshift_dotproduct(g, v_disc, g_obs, v_obs, gp)
    end

    return PointFunction(_internal_turbulent_redshift)
end

# NOTE: this function is not valid within the ISCO, so we need to make sure that
# we **always** set the inner radius of the disc to the ISCO
function velocity_wrapper(m, r, phis, a, M, lum, correlation_length, mach)
    return turbulence_fbmfractal(m, r, phis, a, M, lum, correlation_length, mach)
end

# +------------------------------------------------------------------------------+


# --- Zero Turbulence Thin Disc Line Profile ---
function calculate_line_profile(
    m,
    x,
    d,
    a,
    M,
    lum,
    bins,
    ε,
    turbulenceOn = false,
    correlation_length=1,
    mach=1
    )
    

    if turbulenceOn == false

        _, f = lineprofile(
            bins,
            ε,
            m,
            x,
            d,
            method = BinningMethod(),
            callback = domain_upper_hemisphere(),
            verbose = true
        )

    elseif turbulenceOn == true

        redshift_pf = turbulent_redshift(m, x, velocity_wrapper, a, M, lum, correlation_length, mach)
        pf = redshift_pf ∘ ConstPointFunctions.filter_intersected()
        plane = PolarPlane(GeometricGrid(); Nr = 1000, Nθ = 1000, r_max = outer_radius, r_min = inner_radius)
        _, f = lineprofile(
            ε,
            m,
            x,
            d,
            redshift_pf = pf,
            method = BinningMethod(),
            verbose = true,
            bins = bins,
            plane = plane)

    end

    f[end] = 0
    f[1] = 0

    for i in 1:length(f)
        if f[i] < 0
            f[i] = -f[i]
        end
    end

    return f

end


m = KerrMetric(1.0, 0.998)
inner_radius = Gradus.isco(m)
outer_radius = 400.0

# +--------- Thin Discs ---------+

# Arbitrary Thin Disc
d = ThinDisc(inner_radius, outer_radius)

# +------------------------------+


# +--------- Thick Discs ---------+

# Shakura Sunyaev Disc
#d = ShakuraSunyaev(m, eddington_ratio = 0.3)

# Arbitrary Thick Disc

h_scale = 0.3

function height_profile_opening_angle(ρ, inner_radius, outer_radius, h_scale, angle)
    if ρ < inner_radius || ρ > outer_radius
        return -1.0  # Outside the disc
    else
        p = atan(angle)
        return h_scale * p*ρ
    end
end

function height_profile_pariev_bromley(ρ, inner_radius, outer_radius, h_scale)
    if ρ < inner_radius || ρ > outer_radius
        return -1.0  # Outside the disc
    else
        r_norm = (ρ - inner_radius) / (outer_radius - inner_radius)
        h_factor = h_scale * exp(-r_norm^2)
        return h_factor * ρ
    end
end

"""
d = ThickDisc() do ρ
    height_profile_opening_angle(ρ, inner_radius, outer_radius, h_scale, pi/(0.9))
end
"""

# +------------------------------+

bins = collect(range(0.1, 1.5, 200))
x = SVector(0.0, 1000.0, deg2rad(60), 0.0)
a = 0.998
M = 1.0
lum = 1e46
turbulence_on = true
corona = LampPostModel(h = 10.0)
ε(r) = r^(-7)
f = calculate_line_profile(m, x, d, a, M, lum, bins, ε, turbulence_on, 1)



plot(
    bins, f
)