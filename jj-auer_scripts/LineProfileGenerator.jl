# Script to plot both the zero-turbulence thin disc and perlin/fBm noise turbulence line profiles for different emissivities 
# and inclination angles (emissivities and inclinations angle ranges as of Pariev & Bromley 1998)

# Import libraries
using Plots, Gradus, Measures, LaTeXStrings
include("TurbulenceMaps.jl")
include("SoundSpeed.jl")

plot_font = "Computer Modern"
default(fontfamily=plot_font,
        linewidth=0.8, framestyle=:box, label=nothing, grid=true, yminorgrid = true, xminorgrid = true)
scalefontsizes(1)

incl_angle = 40
x = SVector(0.0, 1_000.0, deg2rad(incl_angle), 0.0)


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
    return turbulence_fbm(m, r, phis, a, M, lum, correlation_length, mach)
end


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

    plane = PolarPlane(GeometricGrid(); Nr = 1000, Nθ = 1000, r_max = 100.0)

    if turbulenceOn == false

        if typeof(ε) <: Gradus.RadialDiscProfile

            _, f = lineprofile(
                m,
                x,
                d,
                ε;
                method = BinningMethod(),
                plane = plane,
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
                method = BinningMethod(),
                plane = plane,
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
                plane = plane,
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
                plane = plane,
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

# Disc geometries
begin
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
end

# How to plot a generic line profile
"""
begin

    # Defaults Parameters:
    begin
        Metric: Kerr
        Disc Type: Shakura-Sunyaev, Eddington Ratio: 0.3
        Radii Limits: ISCO, 100M
        Bin Limits: 0.1, 1.5
        Svector: (0.0, 1000.0, deg2rad(75), 0.0)
        M=1
        Luminosity: 1.9e44 - check if this is 2-10 kev luniosity or bolometric, including all of it 
        Turbulence: Off
        Emissivity: Lamp Post, height = 10M, it's what other people use
        a=0.998
    end

    # Set up figure
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

    # Define parameters
    begin
        a = 0.998
        M = 1.0
        lum = 7.2e42
        lum_edd = LumEdd(M)
        m = KerrMetric(M, a)
        inner_radius = Gradus.isco(m)
        outer_radius = 100.0
        d = ThinDisc(inner_radius, outer_radius)
        bins = collect(range(0.1, 1.5, 200))
        x = SVector(0.0, 1000.0, deg2rad(incl_angle), 0.0)
        turbulenceOn = true
        mach = 1000
        correlation_length=1
    end

    # Calculate flux for each bin and add to plot
    f = calculate_line_profile(m, x, d, a, M, lum, bins, r -> r^-3, turbulenceOn, correlation_length, mach)
    plot!(plt, bins, f, label="Outer Radius = 100M")
    display(plt)
end
"""

# Plotting line profiles for variation of a parameter, for different Mach numbers

fluxes = [
    [[[], [], [], [], []], [[], [], [], [], []], [[], [], [], [], []], [[], [], [], [], []], [[], [], [], [], []]], #Spin
    [[[], [], [], [], []], [[], [], [], [], []], [[], [], [], [], []], [[], [], [], [], []], [[], [], [], [], []]], #Eddington Ratio
    [[[], [], [], [], []], [[], [], [], [], []], [[], [], [], [], []], [[], [], [], [], []], [[], [], [], [], []]], #Correlation Length
    [[[], [], [], [], [], [], []], [[], [], [], [], [], [], []], [[], [], [], [], [], [], []], [[], [], [], [], [], [], []], [[], [], [], [], [], [], []]], #Lamp Height
]

bins = collect(range(0.1, 1.5, 200))

for i in 1:4

    for (j, mach) in enumerate([0, 1, 10, 100, 1000])

        if i == 1 # Spin

            if j == 1

                for (k, a) in enumerate([0.1, 0.5, 0.9, 0.99, 0.998])

                    M = 1.0
                    lum = 7.2e42
                    lum_edd = LumEdd(M)
                    m = KerrMetric(M, a)
                    inner_radius = Gradus.isco(m)
                    outer_radius = 100.0
                    d = ShakuraSunyaev(m, eddington_ratio=0.1)
                    x = SVector(0.0, 1000.0, deg2rad(incl_angle), 0.0)
                    turbulenceOn = false
                    correlation_length=1
                    model = LampPostModel(h=10.0)
                    profile = emissivity_profile(m, d, model)
                    fluxes[i][j][k] = calculate_line_profile(m, x, d, a, M, lum, bins, profile, turbulenceOn, correlation_length, mach)
    
                end
            end

            for (k, a) in enumerate([0.1, 0.5, 0.9, 0.99, 0.998])

                M = 1.0
                lum = 7.2e42
                lum_edd = LumEdd(M)
                m = KerrMetric(M, a)
                inner_radius = Gradus.isco(m)
                outer_radius = 100.0
                d = ShakuraSunyaev(m, eddington_ratio=0.1)
                x = SVector(0.0, 1000.0, deg2rad(incl_angle), 0.0)
                turbulenceOn = true
                correlation_length=1
                model = LampPostModel(h=10.0)
                profile = emissivity_profile(m, d, model)
                fluxes[i][j][k] = calculate_line_profile(m, x, d, a, M, lum, bins, profile, turbulenceOn, correlation_length, mach)

            end
        end

        if i == 2 # Eddington Ratio

            if j == 1

                for (k, ratio) in enumerate([0.1, 0.4, 0.7, 1.0, 1.3])

                    a=0.998
                    M = 1.0
                    lum = 7.2e42
                    lum_edd = LumEdd(M)
                    m = KerrMetric(M, a)
                    inner_radius = Gradus.isco(m)
                    outer_radius = 100.0
                    d = ShakuraSunyaev(m, eddington_ratio=ratio)
                    x = SVector(0.0, 1000.0, deg2rad(incl_angle), 0.0)
                    turbulenceOn = false
                    correlation_length=1
                    model = LampPostModel(h=10.0)
                    profile = emissivity_profile(m, d, model)
                    fluxes[i][j][k] = calculate_line_profile(m, x, d, a, M, lum, bins, profile, turbulenceOn, correlation_length, mach)
    
                end
            end


            for (k, ratio) in enumerate([0.1, 0.4, 0.7, 1.0, 1.3])

                a=0.998
                M = 1.0
                lum = 7.2e42
                lum_edd = LumEdd(M)
                m = KerrMetric(M, a)
                inner_radius = Gradus.isco(m)
                outer_radius = 100.0
                d = ShakuraSunyaev(m, eddington_ratio=ratio)
                x = SVector(0.0, 1000.0, deg2rad(incl_angle), 0.0)
                turbulenceOn = true
                correlation_length=1
                model = LampPostModel(h=10.0)
                profile = emissivity_profile(m, d, model)
                fluxes[i][j][k] = calculate_line_profile(m, x, d, a, M, lum, bins, profile, turbulenceOn, correlation_length, mach)

            end
        end

        if i == 3 # Correlation Length

            if j == 1

                for (k, correlation_length) in enumerate([0.1, 1.0, 5.0, 10.0, 50.0])

                    a=0.998
                    M = 1.0
                    lum = 7.2e42
                    lum_edd = LumEdd(M)
                    m = KerrMetric(M, a)
                    inner_radius = Gradus.isco(m)
                    outer_radius = 100.0
                    d = ShakuraSunyaev(m, eddington_ratio=0.1)
                    x = SVector(0.0, 1000.0, deg2rad(incl_angle), 0.0)
                    turbulenceOn = false
                    model = LampPostModel(h=10.0)
                    profile = emissivity_profile(m, d, model)
                    fluxes[i][j][k] = calculate_line_profile(m, x, d, a, M, lum, bins, profile, turbulenceOn, correlation_length, mach)
    
                end
            end

            for (k, correlation_length) in enumerate([0.1, 1.0, 5.0, 10.0, 50.0])

                a=0.998
                M = 1.0
                lum = 7.2e42
                lum_edd = LumEdd(M)
                m = KerrMetric(M, a)
                inner_radius = Gradus.isco(m)
                outer_radius = 100.0
                d = ShakuraSunyaev(m, eddington_ratio=0.1)
                x = SVector(0.0, 1000.0, deg2rad(incl_angle), 0.0)
                turbulenceOn = true
                model = LampPostModel(h=10.0)
                profile = emissivity_profile(m, d, model)
                fluxes[i][j][k] = calculate_line_profile(m, x, d, a, M, lum, bins, profile, turbulenceOn, correlation_length, mach)

            end
        end

        if i == 4 # Lamp Height

            if j == 1

                for (k, h) in enumerate([2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0])

                    a=0.998
                    M = 1.0
                    lum = 7.2e42
                    lum_edd = LumEdd(M)
                    m = KerrMetric(M, a)
                    inner_radius = Gradus.isco(m)
                    outer_radius = 100.0
                    d = ShakuraSunyaev(m, eddington_ratio=0.1)
                    x = SVector(0.0, 1000.0, deg2rad(incl_angle), 0.0)
                    turbulenceOn = false
                    correlation_length=1
                    model = LampPostModel(h=h)
                    profile = emissivity_profile(m, d, model)
                    fluxes[i][j][k] = calculate_line_profile(m, x, d, a, M, lum, bins, profile, turbulenceOn, correlation_length, mach)
    
                end
            end

            for (k, h) in enumerate([2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0])

                a=0.998
                M = 1.0
                lum = 7.2e42
                lum_edd = LumEdd(M)
                m = KerrMetric(M, a)
                inner_radius = Gradus.isco(m)
                outer_radius = 100.0
                d = ShakuraSunyaev(m, eddington_ratio=0.1)
                x = SVector(0.0, 1000.0, deg2rad(incl_angle), 0.0)
                turbulenceOn = true
                correlation_length=1
                model = LampPostModel(h=h)
                profile = emissivity_profile(m, d, model)
                fluxes[i][j][k] = calculate_line_profile(m, x, d, a, M, lum, bins, profile, turbulenceOn, correlation_length, mach)

            end
        end
    end
end

# Displaying the results
begin
    luxes = fluxes[1][1]
    pal = :seaborn_colorblind

    plt = plot(
        xlabel = L"ν/ν_e \ / \ \textrm{Unitless}",
        ylabel = L"\textrm{Flux \ / \ Arbitrary Units}",
        legend = :topleft,
        left_margin = [5mm 0mm],
        right_margin = [5mm 0mm],
        top_margin = [5mm 0mm],
        bottom_margin = [5mm 0mm],
        palette = pal
        )

    annotate!((1.42, 0.020, (L"\textbf{Laminar}", 10, :black, :center)))
    annotate!((1.42, 0.0185, (L"\mathbf{40 \degree}", 10, :black, :center)))
    
    plot!(plt, bins, luxes[1], label=L"a = 0.1")
    plot!(plt, bins, luxes[2], label=L"a = 0.5")
    plot!(plt, bins, luxes[3], label=L"a = 0.9")
    plot!(plt, bins, luxes[4], label=L"a = 0.99")
    plot!(plt, bins, luxes[5], label=L"a = 0.998")

    display(plt)
end

# Saving current plot
#savefig(plt, "Other/Figs/LeverTweakingPlanned/40/Height/height_1000_40.pdf")