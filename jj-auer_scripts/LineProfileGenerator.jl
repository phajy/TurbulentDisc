# Script to plot both the zero-turbulence thin disc and perlin/fBm noise turbulence line profiles for different emissivities 
# and inclination angles (emissivities and inclinations angle ranges as of Pariev & Bromley 1998)

# Import libraries
using Plots, Gradus, Measures
include("TurbulenceMaps.jl")
include("SoundSpeed.jl")
x = SVector(0.0, 1_000.0, deg2rad(40), 0.0)

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

        if typeof(ε) <: Gradus.RadialDiscProfile

            _, f = lineprofile(
                m,
                x,
                d,
                ε;
                bins = bins,
                verbose = true
            )

        else

            _, f = lineprofile(
                bins,
                ε,
                m,
                x,
                d;
                method = TransferFunctionMethod(),
                verbose = true
            )

        end


    elseif turbulenceOn == true

        redshift_pf = turbulent_redshift(m, x, velocity_wrapper, a, M, lum, correlation_length, mach)
        pf = redshift_pf ∘ ConstPointFunctions.filter_intersected()

        if typeof(ε) <: Gradus.RadialDiscProfile

            _, f = lineprofile(
                m,
                x,
                d,
                ε;
                method = BinningMethod(),
                redshift_pf = pf,
                bins = bins,
                verbose = true)

        else

            _, f = lineprofile(
                bins,
                ε,
                m,
                x,
                d;
                method = BinningMethod(),
                redshift_pf = pf,
                verbose = true)

        end

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




# +--------- Thin Discs ---------+

# Arbitrary Thin Disc
#d = ThinDisc(inner_radius, outer_radius)

# +------------------------------+


# +--------- Thick Discs ---------+

# Shakura Sunyaev Disc


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


"""
Defaults Parameters:

Metric: Kerr
Disc Type: Shakura-Sunyaev, Eddington Ratio: 0.3
Radii Limits: ISCO, 15M
Bin Limits: 0.1, 1.5
Svector: (0.0, 1000.0, deg2rad(40), 0.0)
M=1
Luminosity: 7.2e42 - check if this is 2-10 kev luniosity or bolometric, including all of it 
Turbulence: off
Emissivity: r^(-q), q=3
a=0.998
"""


plt = plot(
    xlabel = "ν / ν_e",
    ylabel = "Flux (Arbitrary Units)",
    title = "Iron Kα Line Profile",
    legend = :topleft,
    left_margin = [5mm 0mm],
    right_margin = [5mm 0mm],
    top_margin = [5mm 0mm],
    bottom_margin = [5mm 0mm]
)

a = 0.998
M = 1.0
lum = 7.2e42
lum_edd = LumEdd(M)
m = KerrMetric(M, a)
inner_radius = Gradus.isco(m)
outer_radius = 15.0
d = Gradus.ThinDisc(inner_radius, outer_radius)
bins = collect(range(0.1, 1.5, 200))
x = SVector(0.0, 1000.0, deg2rad(40), 0.0)
turbulenceOn = false
mach = 1
correlation_length=1
q=3
ε(r) = r^(-q)

f = calculate_line_profile(m, x, d, a, M, lum, bins, ε, turbulenceOn, correlation_length, mach)
plot!(plt, bins, f, ticks=false, legend=false)

display(plt)
savefig(plt, "Other/Figs/LeverTweaking/BogStandardLP.pdf")