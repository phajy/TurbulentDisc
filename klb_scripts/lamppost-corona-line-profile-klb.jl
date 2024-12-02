# Script to plot both the zero-turbulence Shakura-Sunyaev disc and Perlin/fBm noise turbulence line profiles
# using the lamppost corona model

# Import libraries
using Plots, Gradus

# --- Turbulent Shakura-Sunyaev Disc Line Profile ---
include("velocity-structures-klb.jl")

# Parameters for both models
m = KerrMetric(1.0, 0.998)
inner_radius = Gradus.isco(m)
outer_radius = 400.0
# Create Shakura-Sunyaev disc with default parameters
d = ShakuraSunyaev(m)

corona = LampPostModel(h = 10.0)
ε(r) = emissivity_profile(m, d, corona)  # Define emissivity profile using lamppost model

bins = collect(range(0.1, 1.5, 200))
correlation_length = 0.5  # Initial correlation length for turbulence

function turbulent_redshift(metric, x_obs, vel_func, correlation_length)
    g_obs = Gradus.metric(metric, x_obs)
    v_obs = SVector{4, eltype(x_obs)}(1, 0, 0, 0)

    function _internal_turbulent_redshift(m::AbstractMetric, gp, t)
        v_disc = vel_func(m, gp.x[2], gp.x[4], correlation_length)
        g = Gradus.metric(m, gp.x)
        Gradus.RedshiftFunctions._redshift_dotproduct(g, v_disc, g_obs, v_obs, gp)
    end

    return PointFunction(_internal_turbulent_redshift)
end

function velocity_wrapper(m, r, theta, correlation_length)
    return turb_fbm(m, r, theta, correlation_length)
end

function calculate_turbulent_line_profile(bins, m, x, d, correlation_length)
    redshift_pf = turbulent_redshift(m, x, velocity_wrapper, correlation_length)
    pf = redshift_pf ∘ ConstPointFunctions.filter_intersected()
    plane = PolarPlane(GeometricGrid(); Nr = 1000, Nθ = 1000, r_max = outer_radius)
    _, f = lineprofile(bins, m, x, d, redshift_pf = pf, method = BinningMethod(), plane = plane)
    f[end] = 0
    return f
end


# --- Zero Turbulence Shakura-Sunyaev Disc Line Profile ---
        function calculate_zero_turbulence_line_profile(bins, ε, m, x, d)

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
    return f
end


# Define inclination angles
inc_angles = [30, 60, 75]

# Plot for various inclination angles
for inc_angle in inc_angles
    x = SVector(0.0, 1000.0, deg2rad(inc_angle), 0.0)

    # Calculate the zero turbulence line profile
    flux_zero = calculate_zero_turbulence_line_profile(bins, ε, m, x, d)

    # Calculate the turbulent line profile
    flux_turbulent = calculate_turbulent_line_profile(bins, m, x, d, correlation_length)

    # Plot both profiles on the same axes
    plot(
        bins, flux_zero,
        label = "Zero turbulence: i = $inc_angle",
        xlabel = "Redshift",
        ylabel = "Flux (arbitrary units)",
        title = "Shakura-Sunyaev Disc Line Profile",
        legend = :topleft,
        lw = 1
    )
    plot!(
        bins, flux_turbulent,
        linestyle = :dash,
        label = "Turbulent (fBm): i = $inc_angle, L_corr = $correlation_length",
        lw = 1
    )

    # Display the plot for this (q, i) combination
    display(current())
end

