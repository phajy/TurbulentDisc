# Script plots the turbulent velocity using Pariev & Bromley 1998 prescription, as a function of radius in the disc.

# Import relevant libraries and modules
using Plots
using Gradus

# Include the Pariev & Bromley functions script
include("pariev-bromley-equations-klb.jl")

# Define constants
# epsilon = 0.1  # Disc efficiency parameter (example value)
luminosity = 1e46  # Disc luminosity in erg/s (example value)
M_black_hole = 1e8 # Typical mass of black hole (10^8 solar masses)
edd_luminosity_val = eddington_luminosity(M_black_hole) # Eddington luminosity for a 10^8 solar mass black hole
# r_ms = 1.237*M_black_hole # ISCO radius

# Define a range of radii (r) for the plot
radii = collect(range(1.0, stop=25.0, length=200))  # example radial range

# Define values for the spin parameter a/M
spin_parameters = [0.0, 0.5, 0.9, 0.99, 0.998]
m = [KerrMetric(M = 1.0, a = a) for a in spin_parameters]
r_ms = [Gradus.isco(metric) for metric in m]
ϵ = [1.0 - Gradus.CircularOrbits.energy(m[i], r_ms[i]) for i in 1:length(m)]

# Create an empty plot object
plt = plot(
    xlabel = "Radius (r/M)",
    ylabel = "Turbulent Velocity (cₛ / c)",
    title = "Turbulent Velocity vs. Disc Radius",
    legend = :topright,
)

# Calculate and plot turbulent velocity for each spin parameter
for (a, r_ms_val, ϵ_val) in zip(spin_parameters, r_ms, ϵ)
    turbulent_velocities = [sound_speed_ratio(r, a, ϵ_val, luminosity, edd_luminosity_val, r_ms_val) for r in radii]
    plot!(plt, radii, turbulent_velocities, label="a/M = $a")
end

# Show plot
display(plt)
