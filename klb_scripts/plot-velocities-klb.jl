# Import the velocity functions
include("velocity-structures-klb.jl")

# Metric definition
m = KerrMetric(1.0, 0.998)

# Define radii and azimuthal angles
radii = logrange(Gradus.isco(m), 1000, 200)
θ = collect(range(0, 2π, 200))

# Keplerian velocities
keplerian_velocities = Gradus.CircularOrbits.fourvelocity.(m, radii)

# Choose turbulence model by specifying `type` (:rand, :perlin, :rand)
turbulent_velocities = [turbulent_structure(m, r, θ[1]; type=:rand, correlation_length=1) for r in radii]

# --------------------------------------------------------------------- #

# Plot 1: Keplerian vs Turbulent Velocities
begin
    plot(
        radii,
        [i[4] for i in keplerian_velocities],
        xlabel = "r",
        ylabel = "vᵩ",
        xscale = :log10,
        yscale = :log10,
        label = "keplerian",
    )
    plot!(radii, [i[4] for i in turbulent_velocities], label = "turbulent")
end

# Plot 2: Keplerian field in polar projection
begin
    keplerian_field = [
        log10(i[4]) for (i, r) in zip(keplerian_velocities, radii), angle in θ
    ]
    heatmap(θ, radii, keplerian_field, projection=:polar, title="keplerian")
end

# Plot 3: Turbulent field in polar projection
begin
    turbulent_field = [
        log10(turbulent_structure(m, r, angle; type=:rand, correlation_length=10.0)[4]) for r in radii, angle in θ
    ]
    heatmap(θ, radii, turbulent_field, projection=:polar, title="turbulent (rand)")
end

# Plot 4: Comparison map
begin
    comparison_map = [
        (turbulent_structure(m, r, angle; type=:rand, correlation_length=10.0)[4] - v[4]) / v[4]
        for (v, r) in zip(keplerian_velocities, radii), angle in θ
    ]
    heatmap(
        θ,
        radii,
        comparison_map,
        projection=:polar,
        title="(Turbulent - Keplerian) / Keplerian"
    )
end
