# Script to plot thin disc line profiles for different emissivities and inclination angles
# (refer to Fig 5 Pariev & Bromley 1998)

using Plots, Gradus

# Define inclination angles and emissivity indices
inc_angles = [30, 60, 75]
q_values = [2, 3, 4]

# Create the thin disc with inner and outer radius
inner_radius = 0.0
outer_radius = 400.0
d = ThinDisc(inner_radius, outer_radius)

# Kerr metric setup
m = KerrMetric(1.0, 0.998)

# Define a custom bin range for redshift (g grid for line profile)
bins = collect(range(0.1, 1.5, 200))

# Function to calculate line profile for the thin disc
function calculate_line_profile(m, x, d, bins, q)
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

# Loop through each combination of emissivity index and inclination angle and plot
for q in q_values
    for inc_angle in inc_angles
        # Observer position vector with inclination angle
        x = SVector(0.0, 1000.0, deg2rad(inc_angle), 0.0)

        # Calculate the line profile flux for the current parameters
        flux = calculate_line_profile(m, x, d, bins, q)

        # Plot the line profile with annotations for inclination angle and emissivity index
        plot(
            bins, flux,
            xlabel = "Redshift",
            ylabel = "Flux (arbitrary units)",
            legend = false,
            title = "i = $inc_angle, q = $q"
        )
        
        # Display each plot
        display(current())
    end
end
