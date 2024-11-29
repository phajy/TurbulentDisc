# Import the velocity functions
include("velocity-structures-klb.jl")

# Metric definition
m = KerrMetric(1.0, 0.998)

# Define radii and azimuthal angles
radii = logrange(Gradus.isco(m), 1000, 200)
θ = collect(range(0, 2π, length=200))

# Keplerian velocities
keplerian_velocities = Gradus.CircularOrbits.fourvelocity.(m, radii)

# Select turbulence type (choose one: :random, :perlin, :fbm)
selected_type = :random

# Function to generate turbulent field data for plotting
function generate_field(m, radii, θ, type, correlation_length)
    return [
        log10(turbulent_structure(m, r, θ_val; type=type, correlation_length=correlation_length)[4])
        for r in radii, θ_val in θ
    ]
end

# Function to generate a single comparison map
function generate_comparison_map(m, radii, θ, keplerian_velocities, type, correlation_length)
    return [
        (turbulent_structure(m, r, θ_val; type=type, correlation_length=correlation_length)[4] - v[4]) / v[4]
        for (v, r) in zip(keplerian_velocities, radii), θ_val in θ
    ]
end

# --------------------------------------------------------------------- #

# Generate Keplerian Field
keplerian_field = [
    log10(i[4]) for (i, r) in zip(keplerian_velocities, radii), θ_val in θ
]

# Generate Turbulent Field for the Selected Type
turbulent_field = generate_field(m, radii, θ, selected_type, 1)

# Generate Comparison Map
comparison_map = generate_comparison_map(m, radii, θ, keplerian_velocities, selected_type, 1)

# --------------------------------------------------------------------- #

# Plot Keplerian Field
heatmap(
    θ,
    radii,
    keplerian_field,
    projection=:polar,
    title="Keplerian Field",
)

# Plot Turbulent Field
heatmap(
    θ,
    radii,
    turbulent_field,
    projection=:polar,
    title="Turbulent Field ($selected_type)",
)

# Plot Comparison Map
heatmap(
    θ,
    radii,
    comparison_map,
    projection=:polar,
    title="Comparison: (Turbulent - Keplerian) / Keplerian ($selected_type)",
)
