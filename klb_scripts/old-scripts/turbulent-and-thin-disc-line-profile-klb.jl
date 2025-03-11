# Script to plot both the zero-turbulence thin disc and Perlin/fBm noise turbulence line profiles 
# for different emissivities and inclination angles.

# Import libraries
using Plots, Gradus, LaTeXStrings

include("velocity-structures-klb.jl")
include("pariev-bromley-equations-klb.jl")

# Choose turbulence model: "perlin" or "fbm" 
turbulence_model = "fbm" # Options: "perlin" or "fbm"

# Choose emissivity model: "powerlaw" or "lamppost"
emissivity_model = "powerlaw"  # Options: "powerlaw" or "lamppost" 

# --- Turbulent Thin Disc Line Profile ---
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

function velocity_wrapper(m, r, theta, correlation_length, a, M, L, L_edd, r_ms, epsilon, mach)
    if turbulence_model == "perlin"
        return turb_perlin(m, r, theta, a, M, L, L_edd, r_ms, epsilon, correlation_length, mach)
    elseif turbulence_model == "fbm"
        return turb_fbm(m, r, theta, a, M, L, L_edd, r_ms, epsilon, correlation_length, mach)
    end
end

function get_emissivity_function(m, d, q, height, emissivity_model)
    if emissivity_model == "powerlaw"
        return r -> r^(-q)
    elseif emissivity_model == "lamppost"
        corona = LampPostModel(h = height)
        return emissivity_profile(m, d, corona)
    else
        error("Invalid emissivity model. Choose 'powerlaw' or 'lamppost'.")
    end
end

function calculate_turbulent_line_profile(m, x, d, bins, correlation_length, q, a, M, L, L_edd, r_ms, epsilon, mach, emissivity_model)
    redshift_pf = turbulent_redshift(m, x, velocity_wrapper, correlation_length, a, M, L, L_edd, r_ms, epsilon, mach)
    pf = redshift_pf ∘ ConstPointFunctions.filter_intersected()
    plane = PolarPlane(GeometricGrid(); Nr = 1000, Nθ = 1000, r_max = outer_radius)
    ε = get_emissivity_function(m, d, q, height, emissivity_model)
    if emissivity_model == "powerlaw"
        _, f = lineprofile(bins, ε, m, x, d, redshift_pf = pf, method = BinningMethod(), plane = plane)
    elseif emissivity_model == "lamppost"
        _, f = lineprofile(m, x, d, ε; bins, redshift_pf = pf, method = BinningMethod(), plane = plane)
    end
    f[end] = 0
    return f
end

# --- Zero Turbulence Thin Disc Line Profile ---
function calculate_zero_turbulence_line_profile(m, x, d, bins, q, emissivity_model)
    ε = get_emissivity_function(m, d, q, height, emissivity_model)
    _, f = lineprofile(bins, ε, m, x, d, method = BinningMethod(), callback = domain_upper_hemisphere(), verbose = true)
    return f
end

# Parameters
m = KerrMetric(1.0, 0.998)
inner_radius = Gradus.isco(m)
outer_radius = 15.0
bins = collect(range(0.1, 1.5, 200))
correlation_length = 1.0
L_eddington(M) = 1.2e46 * (M/1e8)

# Define inclination angles, emissivity indices, and other parameters
inc_angles = [30, 60, 75]
q_values = [2, 3, 4]
a = 0.998
M = 1.0
L = 1e46
L_edd = L_eddington(M)
r_ms = Gradus.isco(m)
epsilon = 0.1
mach = 1
height = 10.0  # Height for lamppost model

d = ThinDisc(inner_radius, outer_radius)  # Define thin disc

# Loop over parameters and generate plots
for q in q_values
    for inc_angle in inc_angles
        x = SVector(0.0, 1000.0, deg2rad(inc_angle), 0.0)

        # Calculate zero turbulence profile
        flux_zero = calculate_zero_turbulence_line_profile(m, x, d, bins, q, emissivity_model)

        # Calculate turbulent profile
        flux_turbulent = calculate_turbulent_line_profile(m, x, d, bins, correlation_length, q, a, M, L, L_edd, r_ms, epsilon, mach, emissivity_model)

        plot(
            bins, flux_zero,
            label = "Zero turbulence",
            xlabel = "Redshift",
            ylabel = "Flux (arbitrary units)",
            legend = :topleft,
            lw = 0.8,
            color = :black,
            linestyle = :solid
        )
        plot!(
            bins, flux_turbulent,
            label = "Turbulent " * L"(M = %$mach)",
            lw = 1, 
            color = :black,
            linestyle = :dashdot
        )

        annot_text_i = L"i = %$inc_angle^{\circ}"
        annot_text_q = L"q = %$q"

        annot_x = 1.4
        annot_y_max = maximum(flux_zero)
        offset = 0.1 * annot_y_max 
        
        annot_y_i = annot_y_max - offset 
        annot_y_q = annot_y_i - (1.1*offset)    

        annotate!(annot_x, annot_y_i, text(annot_text_i, 11, :black, :right))
        annotate!(annot_x, annot_y_q, text(annot_text_q, 11, :black, :right))

        display(current())
        
        # Save plot
        output_dir = "C:\\Users\\Kate\\project\\TurbulentDisc\\klb_plots\\line-profiles-thin\\$(turbulence_model)-line-profiles\\$(emissivity_model)"
        mkpath(output_dir)
        filename = joinpath(output_dir, "line_profile_$(emissivity_model)_i$(inc_angle)_q$(q)_M$(mach).png")
        savefig(filename)
        println("Saved figure to: $filename")
    end 
end
