# Script to plot both the zero-turbulence Shakura-Sunyaev disc and Perlin/fBm noise turbulence line profiles 
# for different emissivities and inclination angles.

# Import libraries
using Plots, Gradus, LaTeXStrings

include("velocity-structures-klb.jl")
include("pariev-bromley-equations-klb.jl")

L_eddington(M) = 1.2e46 * (M/1e8) # Eddington luminosity

# --- Turbulent Shakura-Sunyaev Disc Line Profile ---
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
    return turb_perlin(m, r, theta, a, M, L, L_edd, r_ms, epsilon, correlation_length, mach)
end

# Function to calculate emissivity based on chosen model
function get_emissivity_function(m, d, q, emissivity_model)
    if emissivity_model == "powerlaw"
        return r -> r^(-q)
    elseif emissivity_model == "lamppost"
        corona = LampPostModel(h = 10.0)  # Define lamp-post corona
        return emissivity_profile(m, d, corona)
    else
        error("Invalid emissivity model. Choose 'powerlaw' or 'lamppost'.")
    end
end

function calculate_turbulent_line_profile(m, x, d, bins, correlation_length, q, a, M, L, L_edd, r_ms, epsilon, mach, emissivity_model)
    redshift_pf = turbulent_redshift(m, x, velocity_wrapper, correlation_length, a, M, L, L_edd, r_ms, epsilon, mach)
    pf = redshift_pf ∘ ConstPointFunctions.filter_intersected()
    plane = PolarPlane(GeometricGrid(); Nr = 1000, Nθ = 1000, r_max = outer_radius)
    ε = get_emissivity_function(m, d, q, emissivity_model)
    _, f = lineprofile(bins, ε, m, x, d, redshift_pf = pf, method = BinningMethod(), plane = plane)
    f[end] = 0
    return f
end

# --- Zero Turbulence Shakura-Sunyaev Disc Line Profile ---
function calculate_zero_turbulence_line_profile(m, x, d, bins, q, emissivity_model)
    ε = get_emissivity_function(m, d, q, emissivity_model)
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

# --- Combined Plotting ---
# Parameters for both models
m = KerrMetric(1.0, 0.998)
inner_radius = Gradus.isco(m)
outer_radius = 400.0
bins = collect(range(0.1, 1.5, 200))
correlation_length = 1  # Correlation length for turbulence

# Define inclination angles and emissivity indices
inc_angles = [30, 60, 75]
q_values = [2, 3, 4]

# Choose emissivity model: "powerlaw" or "lamppost"
emissivity_model = "powerlaw"  # or "lamppost" 

# Define parameters for sound speed normalisation
a = 0.998  # Black hole spin
M = 1.0  # Black hole mass in geometrized units
L = 1e46  # Luminosity in ergs/s
L_edd = L_eddington(M)  # Compute Eddington luminosity
r_ms = Gradus.isco(m)  # Compute ISCO
epsilon = 0.1  # Efficiency factor
mach = 1 # Mach number

# Plot for each combination of parameters
for q in q_values
    for inc_angle in inc_angles
        x = SVector(0.0, 1000.0, deg2rad(inc_angle), 0.0)

        # Calculate the zero turbulence line profile (set Eddington ratio to 0.0 for a truly thin disc)
        flux_zero = calculate_zero_turbulence_line_profile(m, x, ShakuraSunyaev(m, eddington_ratio=0.3), bins, q, emissivity_model)

        # Calculate the turbulent line profiles for different Eddington ratios
        flux_turbulent_05 = calculate_turbulent_line_profile(m, x, ShakuraSunyaev(m, eddington_ratio=0.5), bins, correlation_length, q, a, M, L, L_edd, r_ms, epsilon, mach, emissivity_model)
        flux_turbulent_1 = calculate_turbulent_line_profile(m, x, ShakuraSunyaev(m, eddington_ratio=1.0), bins, correlation_length, q, a, M, L, L_edd, r_ms, epsilon, mach, emissivity_model)

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
            bins, flux_turbulent_05,
            label = "Turbulent " * L"(L_{Edd} = 0.5, M = %$mach)",
            lw = 1, 
            color = :black,
            linestyle = :dashdot
        )
        plot!(
            bins, flux_turbulent_1,
            label = "Turbulent " * L"(L_{Edd} = 1.0, M = %$mach)",
            lw = 1, 
            color = :black,
            linestyle = :dash
        )

        annot_text_i = L"i = %$inc_angle^{\circ}"
        if emissivity_model == "powerlaw"
            annot_text_q = L"q = %$q"
        end

        annot_x = 1.4

        annot_y_max = maximum(flux_zero)
        offset = 0.1 * annot_y_max 
        
        annot_y_i = annot_y_max - offset 
        annot_y_q = annot_y_i - (1.1*offset)    

        annotate!(annot_x, annot_y_i, text(annot_text_i, 11, :black, :right))
        if emissivity_model == "powerlaw"
            annotate!(annot_x, annot_y_q, text(annot_text_q, 11, :black, :right))
        end

        # Display the plot for this (q, i) combination
        display(current())
        
        
        # Save plot for this (q, i) combination
        #output_dir = raw"C:\Users\Kate\project\TurbulentDisc\klb_plots\line-profiles\perlin-line-profiles"
        #mkpath(output_dir)

        # Generate filename with inclination angle, emissivity index, and Mach number
        #filename = joinpath(output_dir, "line_profile_i$(inc_angle)_q$(q)_M$(mach).png")

        # Save the figure
        #savefig(filename)
        #println("Saved figure to: $filename")

        

    end
end
