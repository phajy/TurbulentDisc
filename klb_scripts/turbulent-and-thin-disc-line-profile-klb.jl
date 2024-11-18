# Script to plot both the zero-turbulence thin disc and perlin noise turbulence line profiles for different emissivities 
# and inclination angles (emissivities and inclinations angle ranges as of Pariev & Bromley 1998)

# Import libraries
using Plots, Gradus

# --- Turbulent Thin Disc Line Profile ---
include("velocity-structures-klb.jl")

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

function calculate_turbulent_line_profile(m, x, d, bins, correlation_length)
    redshift_pf = turbulent_redshift(m, x, velocity_wrapper, correlation_length)
    pf = redshift_pf ∘ ConstPointFunctions.filter_intersected()
    plane = PolarPlane(GeometricGrid(); Nr = 1000, Nθ = 1000, r_max = 5 * d.outer_radius)
    _, f = lineprofile(m, x, d, redshift_pf = pf, method = BinningMethod(), bins = bins, plane = plane)
    f[end] = 0
    return f
end


# --- Zero Turbulence Thin Disc Line Profile ---
function calculate_zero_turbulence_line_profile(m, x, d, bins, q)
    ε(r) = r^(-q)  # Define emissivity function with given index q
    _, f = lineprofile(
        m,
        x,
        d,
        method = BinningMethod(),
        callback = domain_upper_hemisphere(),
        verbose = true,
        bins = bins
    )
    return f
end


# --- Combined Plotting ---
# Parameters for both models
m = KerrMetric(1.0, 0.998)
inner_radius = Gradus.isco(m)
outer_radius = 400.0
d = ThinDisc(inner_radius, outer_radius)
bins = collect(range(0.1, 1.5, 200))
correlation_length = 1.0  # Initial correlation length for turbulence

# Define inclination angles and emissivity indices
inc_angles = [30, 60, 75]
q_values = [2, 3, 4]

# Plot for each combination of parameters
for q in q_values
    for inc_angle in inc_angles
        x = SVector(0.0, 1000.0, deg2rad(inc_angle), 0.0)

        # Calculate the zero turbulence line profile
        flux_zero = calculate_zero_turbulence_line_profile(m, x, d, bins, q)

        # Calculate the turbulent line profile
        flux_turbulent = calculate_turbulent_line_profile(m, x, d, bins, correlation_length)

        # Plot both profiles on the same axes
        plot(
            bins, flux_zero,
            label = "Zero turbulence: i = $inc_angle, q = $q",
            xlabel = "Redshift",
            ylabel = "Flux (arbitrary units)",
            legend = :topleft,
            lw = 1
        )
        plot!(
            bins, flux_turbulent,
            linestyle = :dash,
            label = "Turbulent: i = $inc_angle, q = $q",
            lw = 1
        )

        # Display the plot for this (q, i) combination
        display(current())
    end
end
