# Script to plot the zero-turbulence and perlin/fBm noise turbulence line profiles for a arbitrarily thick disc for different emissivities 
# and inclination angles (emissivities and inclinations angle ranges as of Pariev & Bromley 1998)

# Import relevant libraries
using Plots, Gradus

# Import turbulent velocity functions
include("velocity-structures-klb.jl")
include("pariev-bromley-equations-klb.jl")

# Define the Kerr metric
m = Gradus.KerrMetric(M=1.0, a=0.998)

# Define the radii
inner_radius = Gradus.isco(m)
outer_radius = 400.0

L_eddington(M) = 1.2e46 * (M/1e8) # edd luminosity

# Define the disc height profile, this is handled the same way as the finite thickness disc in 
# Pariev & Bromley (1998), which is based on Novikov & Thorne (1973) (ultimately the Shakura-Sunyaev model,
# where h_scale < 1)
function thick_disc_height_profile(ρ, inner_radius, outer_radius, h_scale)
    if ρ < inner_radius || ρ > outer_radius
        return -1.0  # Outside the disc
    else
        r_norm = (ρ - inner_radius) / (outer_radius - inner_radius)
        h_factor = h_scale * exp(-r_norm^2)
        return h_factor * ρ
    end
end

# Create the thick disc geometry
h_scale = 0.2
thick_disc = Gradus.ThickDisc() do ρ
    thick_disc_height_profile(ρ, inner_radius, outer_radius, h_scale)
end

# Define inclination angles and emissivity indices
inc_angles = [30, 60, 75]
q_values = [2, 3, 4]

# Define redshift bins and turbulence parameters
bins = collect(range(0.1, 1.5, length=200))
correlation_length = 1.0

# Define parameters for sound speed normalisation
a = 0.998  # Black hole spin
M = 1.0  # Black hole mass in geometrised units
L = 1e46  # Luminosity in ergs/s
L_edd = L_eddington(M)  # Compute Eddington luminosity
r_ms = Gradus.isco(m)  # Compute ISCO
epsilon = 0.1  # Efficiency factor

# Turbulent redshift function
function turbulent_redshift(metric, x_obs, vel_func, correlation_length, a, M, L, L_edd, r_ms, epsilon)
    g_obs = Gradus.metric(metric, x_obs)
    v_obs = SVector{4, eltype(x_obs)}(1, 0, 0, 0)

    function _internal_turbulent_redshift(m::AbstractMetric, gp, t)
        v_disc = vel_func(m, gp.x[2], gp.x[4], correlation_length, a, M, L, L_edd, r_ms, epsilon)
        g = Gradus.metric(m, gp.x)
        Gradus.RedshiftFunctions._redshift_dotproduct(g, v_disc, g_obs, v_obs, gp)
    end

    return PointFunction(_internal_turbulent_redshift)
end

function velocity_wrapper(m, r, theta, correlation_length, a, M, L, L_edd, r_ms, epsilon)
    return turb_fbm(m, r, theta, correlation_length, a, M, L, L_edd, r_ms, epsilon)
end

# Line profile functions
function calculate_turbulent_line_profile(m, x, d, bins, correlation_length, a, M, L, L_edd, r_ms, epsilon)
    redshift_pf = turbulent_redshift(m, x, velocity_wrapper, correlation_length, a, M, L, L_edd, r_ms, epsilon)
    pf = redshift_pf ∘ ConstPointFunctions.filter_intersected()
    plane = PolarPlane(GeometricGrid(); Nr = 1000, Nθ = 1000, r_max = outer_radius, r_min = inner_radius)
    ε(r) = r^(-q)
    _, f = lineprofile(
        bins,      
        ε,   
        m,
        x,
        d,
        redshift_pf = pf,
        method = BinningMethod(),
        plane = plane
    )
    f[end] = 0
    return f
end

function calculate_zero_turbulence_line_profile(m, x, d, bins, q)
    ε(r) = r^(-q)
    _, f = lineprofile(
        bins,      
        ε,   
        m,
        x,
        d,
        method = BinningMethod(),
        callback = domain_upper_hemisphere(),
        verbose = true,
    )
    return f
end

# Iterate through all combinations of inclination angles and emissivity indices
for q in q_values
    for i in inc_angles
        x_obs = SVector(0.0, 1000.0, deg2rad(i), 0.0)
        
        # Compute line profiles
        flux_zero_turbulence = calculate_zero_turbulence_line_profile(m, x_obs, thick_disc, bins, q)
        flux_turbulent = calculate_turbulent_line_profile(m, x_obs, thick_disc, bins, correlation_length, a, M, L, L_edd, r_ms, epsilon)
        
        # Create a new plot for each combination
        plot(
            bins, flux_zero_turbulence,
            label = "Zero turb: i=$i, q=$q",
            lw = 1, color = :blue,
            xlabel = "Redshift", ylabel = "Flux",
            title = "Thick Disc Line Profile",
            legend = :topleft
        
        )
        
        # Add turbulent profile to the same plot
        plot!(
            bins, flux_turbulent,
            label = "Turbulent (fBm): i=$i, q=$q, L_corr = $correlation_length",
            lw = 1, linestyle = :dash, color = :red
        )
        
        # Display the plot for this (q, i) combination
        display(current())
    end
end
