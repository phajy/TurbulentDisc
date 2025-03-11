# Import required libraries
using Plots, Gradus, LaTeXStrings

# Import the velocity functions
include("velocity-structures-klb.jl")
include("pariev-bromley-equations-klb.jl")

# Define correlation length
correlation_length = 0.1
turbulence_model = "perlin"  # "perlin" or "fbm"

# Define inclination angles and emissivity indices
inc_angles = [30, 60, 75]
q_values = [2, 3, 4]

# Function to calculate non-turbulent (Keplerian) redshift
function keplerian_redshift(metric, x_obs)
    g_obs = Gradus.metric(metric, x_obs)  # Metric at observer's position
    v_obs = SVector{4, eltype(x_obs)}(1, 0, 0, 0)  # Stationary observer velocity

    function _internal_keplerian_redshift(m::AbstractMetric, gp, t)
        v_disc = Gradus.CircularOrbits.fourvelocity(m, gp.x[2])  # Keplerian four-velocity
        g = Gradus.metric(m, gp.x)
        Gradus.RedshiftFunctions._redshift_dotproduct(g, v_disc, g_obs, v_obs, gp)
    end

    return PointFunction(_internal_keplerian_redshift)
end

# Function to calculate turbulent redshift
function turbulent_redshift(metric, x_obs, vel_func, correlation_length, a, M, L, L_edd, r_ms, epsilon, mach)
    g_obs = Gradus.metric(metric, x_obs)
    v_obs = SVector{4, eltype(x_obs)}(1, 0, 0, 0)

    function _internal_turbulent_redshift(m::AbstractMetric, gp, t)
        v_disc = vel_func(m, gp.x[2], gp.x[4], correlation_length, a, M, L, L_edd, r_ms, epsilon, mach)
        g = Gradus.metric(m, gp.x)
        Gradus.RedshiftFunctions._redshift_dotproduct(g, v_disc, g_obs, v_obs, gp)
    end

    return PointFunction(_internal_turbulent_redshift)
end

# Wrapper function to select the appropriate turbulence model
function velocity_wrapper(m, r, theta, correlation_length, a, M, L, L_edd, r_ms, epsilon, mach)
    if turbulence_model == "perlin"
        return turb_perlin(m, r, theta, a, M, L, L_edd, r_ms, epsilon, correlation_length, mach)
    elseif turbulence_model == "fbm"
        return turb_fbm(m, r, theta, a, M, L, L_edd, r_ms, epsilon, correlation_length, mach)
    end
end

# Parameters and initial setup
m = KerrMetric(1.0, 0.998)
inner_radius = Gradus.isco(m)
outer_radius = 400.0

# Define parameters for turbulence model
a = 0.998  # Black hole spin
M = 1.0  # Black hole mass in geometrized units
L = 1e46  # Luminosity in ergs/s
L_edd = 1.2e46 * (M/1e8)  # Compute Eddington luminosity
r_ms = Gradus.isco(m)  # Compute ISCO
epsilon = 0.1  # Efficiency factor
mach = 1  # Mach number

# Define the Shakura-Sunyaev disc
eddington_ratios = [0.3, 0.5, 1.0]
discs = Dict(0.3 => ShakuraSunyaev(m, eddington_ratio=0.3),
             0.5 => ShakuraSunyaev(m, eddington_ratio=0.5),
             1.0 => ShakuraSunyaev(m, eddington_ratio=1.0))

# Iterate over inclination angles, emissivity indices, and Eddington ratios
for q in q_values
    for inc_angle in inc_angles
        for edd_ratio in eddington_ratios
            x = SVector(0.0, 1000.0, deg2rad(inc_angle), 0.0)
            d = discs[edd_ratio]

            # Compute redshift functions
            keplerian_pf = keplerian_redshift(m, x) ∘ ConstPointFunctions.filter_intersected()
            redshift_pf = turbulent_redshift(m, x, velocity_wrapper, correlation_length, a, M, L, L_edd, r_ms, epsilon, mach)
            turbulent_pf = redshift_pf ∘ ConstPointFunctions.filter_intersected()

            # Generate non-turbulent (keplerian) redshift image
            α, β, img_keplerian = rendergeodesics(
                m, x, d, 20_000.0,
                αlims=(-12, 12), βlims=(-8, 8),
                image_width=800, image_height=400,
                verbose=true, pf=keplerian_pf
            )

            # Generate turbulent redshift image
            _, _, img_turbulent = rendergeodesics(
                m, x, d, 20_000.0,
                αlims=(-12, 12), βlims=(-8, 8),
                image_width=800, image_height=400,
                verbose=true, pf=turbulent_pf
            )

            # Compute the difference between turbulent and non-turbulent images
            Δz = img_turbulent .- img_keplerian

            # Define output directory
            output_dir = "C:\\Users\\Kate\\project\\TurbulentDisc\\klb_plots\\ray-traced\\$(turbulence_model)-ray-traced\\keplerian-turbulent-difference"
            mkpath(output_dir)

            # Plot non-turbulent (keplerian) image
            heatmap(α, β, img_keplerian, aspect_ratio=1, title="Non-Turbulent Redshift Image")

            # Plot turbulent image
            heatmap(α, β, img_turbulent, aspect_ratio=1, title="Turbulent Redshift Image")

            # Plot and save difference image
            filename_difference = joinpath(output_dir, "redshift_difference_i$(inc_angle)_q$(q)_Ledd$(edd_ratio).png")
            heatmap(α, β, Δz, aspect_ratio=1, title="Redshift Difference (Turbulence Effect)", color=:balance)
            savefig(filename_difference)

        end
    end
end
