# import functions
include("TurbulenceMaps.jl")

m = KerrMetric(1.0, 0.998)

# log range from the innermost stable circular orbit (inner edge of the
# accretion disc) out to some arbitrary radius
radii = logrange(Gradus.isco(m), 1000, 200)

# if you need it, these will be the azimuthal angles to each "radial" bin
θ = collect(range(0, 2π, 200))

keplerian_velocities = Gradus.CircularOrbits.fourvelocity.(m, radii)

a = 0.998
M = 1.0
lum = 7.2e42
correlation_length = 1
mach = 1

# Calculate turbulent velocities for each combination of r and θ
turbulent_velocities = [
    [turbulence_fbmfractal(m, radii[i], θ[j], a, M, lum, correlation_length, mach) for j in 1:length(θ)]
    for i in 1:length(radii)
]

# Calculate the difference between turbulent and keplerian velocities
velocity_differences = [
    [turbulent_velocities[i][j] - keplerian_velocities[i] for j in 1:length(θ)]
    for i in 1:length(radii)
]

# Calculate magnitudes of the differences (excluding the first component)
difference_mags = [
    [sqrt(sum(x -> x^2, velocity_differences[i][j][2:4])) for j in 1:length(θ)]
    for i in 1:length(radii)
]

print(keplerian_velocities[1][4])
print(turbulent_velocities[1][1][4])
turbulence_fbmfractal(m, radii[20], θ[20], a, M, lum, correlation_length, mach)

# --------------------------------------------------------------------- #

# some visualisation methods (use the ones that make the most sense to you)
# note these are all only using the 4th component of the velocity vector (which
# is the azimuthal component):
# 1 - time
# 2 - radial (r)
# 3 - poloidal (θ)
# 4 - azimuthal (ϕ)

"""
========================================================
# This is a 1D plot of azimuthal velocity against radius
========================================================

begin
    plot(
        radii,
        [sqrt(sum(x -> x^2, v)) for v in keplerian_velocities],
        xlabel = "r",
        ylabel = "Velocity Magnitude",
        xscale = :log10,
        yscale = :log10,
        label = "keplerian",
    )
    plot!(radii, [sqrt(sum(x -> x^2, turbulent_velocities[i][1])) for i in 1:length(radii)], label = "turbulent")
end


begin
    keplerian_field =
        [sqrt(sum(x -> x^2, v)) for v in keplerian_velocities, angle in θ]
    heatmap(θ, radii, keplerian_field, projection = :polar, title = "keplerian")
end

begin
    turbulent_field = [sqrt(sum(x -> x^2, turbulent_velocities[i][j])) for i in 1:length(radii), j in 1:length(θ)]
    heatmap(θ, radii, turbulent_field, projection = :polar, title = "turbulent")
end
"""

begin
    difference_field = [difference_mags[i][j] for i in 1:length(radii), j in 1:length(θ)]
    heatmap(
        θ,
        radii,
        difference_field,
        projection = :polar,
        title = "Magnitude of Velocity Difference",
    )
end