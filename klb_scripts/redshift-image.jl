# import required libraries
using Plots, Gradus, LaTeXStrings

# import the velocity functions
include("velocity-structures-klb.jl")
include("pariev-bromley-equations-klb.jl")

# Define correlation length
correlation_length = 10

turbulence_model = "perlin" # "perlin" or "fbm"

# Define inclination angles and emissivity indices
inc_angles = [30, 60, 75]
q_values = [2, 3, 4]

function turbulent_redshift(metric, x_obs, vel_func, correlation_length, a, M, L, L_edd, r_ms, epsilon, mach)
    # Metric matrix at the observer's position
    g_obs = Gradus.metric(metric, x_obs)
    # Fixed stationary observer velocity
    v_obs = SVector{4, eltype(x_obs)}(1, 0, 0, 0)

    # Internal closure function
    function _internal_turbulent_redshift(m::AbstractMetric, gp, t)
        # Pass correlation_length to vel_func
        v_disc = vel_func(m, gp.x[2], gp.x[4], correlation_length, a, M, L, L_edd, r_ms, epsilon, mach)

        g = Gradus.metric(m, gp.x)
        Gradus.RedshiftFunctions._redshift_dotproduct(g, v_disc, g_obs, v_obs, gp)
    end

    return PointFunction(_internal_turbulent_redshift)
end

# Wrapper function to select appropriate turbulence model
function velocity_wrapper(m, r, theta, correlation_length, a, M, L, L_edd, r_ms, epsilon, mach)
    if turbulence_model == "perlin"
        return turb_perlin(m, r, theta, a, M, L, L_edd, r_ms, epsilon, correlation_length, mach)
    elseif turbulence_model == "fbm"
        return return turb_fbm(m, r, theta, a, M, L, L_edd, r_ms, epsilon, correlation_length, mach)
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

# Iterate over inclination angles and emissivity indices
for q in q_values
    for inc_angle in inc_angles
        for edd_ratio in eddington_ratios
            x = SVector(0.0, 1000.0, deg2rad(inc_angle), 0.0)
            d = discs[edd_ratio]

            # Compose our turbulent redshift function with a filter function to remove points outside of the ISCO
            redshift_pf = turbulent_redshift(m, x, velocity_wrapper, correlation_length, a, M, L, L_edd, r_ms, epsilon, mach)
            pf = redshift_pf ∘ ConstPointFunctions.filter_intersected()

            # Render geodesics using the modified point function 
            α, β, img = rendergeodesics(
                m,
                x,
                d,
                20_000.0,  # Maximum integration time
                αlims = (-12, 12), 
                βlims = (-8, 8),
                image_width = 800,
                image_height = 400,
                verbose = true,
                pf = pf,
            )


            # Plot ray trace geodesic image
            heatmap(α, β, img, aspect_ratio = 1, xlabel="α", ylabel="β", title="Redshift Image: i=$(inc_angle)°, q=$(q), L_Edd=$(edd_ratio), corr=$(correlation_length)")

            # Save geodesic image
            #output_dir = "C:\\Users\\Kate\\project\\TurbulentDisc\\klb_plots\\ray-traced\\$(turbulence_model)-ray-traced\\powerlaw"
            output_dir = "C:\\Users\\Kate\\project\\TurbulentDisc\\klb_plots\\ray-traced\\$(turbulence_model)-ray-traced\\powerlaw\\zoomed-in"
            mkpath(output_dir)
            filename = joinpath(output_dir, "ray_traced_powerlaw_i$(inc_angle)_q$(q)_Ledd$(edd_ratio)_M$(mach)_corr$(correlation_length).png")
            savefig(filename)
            println("Saved figure to: $filename")

            # Define emissivity function 
            ϵ(r) = r^(-q)

            # Create and plot line profile 
            bins = collect(range(0.0, 2.0, 200))
            plane = PolarPlane(GeometricGrid(); Nr = 1000, Nθ = 1000, r_max = outer_radius)
            bins, f = lineprofile(bins, ϵ, m, x, d, redshift_pf = pf, verbose = true, method = BinningMethod(), plane = plane)


            # Set whatever is in the last bin to 0 as it's most likely a noise contribution
            f[end] = 0

            # Plot line profile
            plot(bins, f, 
                legend = false,
                xlabel="Redshift",
                ylabel="Flux (arbitrary units)",
                title="Line Profile: i=$(inc_angle)°, q=$(q), L_Edd=$(edd_ratio)")
        end
    end
end