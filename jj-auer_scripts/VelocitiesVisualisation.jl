# import functions
using Gradus, Plots
include("TurbulenceMaps.jl")

m = KerrMetric(1.0, 0.998)

# log range from the innermost stable circular orbit (inner edge of the
# accretion disc) out to some arbitrary radius
radii = logrange(Gradus.isco(m), 5, 200)

# if you need it, these will be the azimuthal angles to each "radial" bin
θ = collect(range(0, 2π, 200))

keplerian_velocities = Gradus.CircularOrbits.fourvelocity.(m, radii)

a = 0.998
M = 1.0
lum = 7.2e42
correlation_length = 1
mach = 1000


function velocity_wrapper(m, r, phis, a, M, lum, correlation_length, mach)
    return turbulence_perlin(m, r, phis, a, M, lum, correlation_length, mach)
end

# Calculate turbulent velocities for each combination of r and θ
turbulent_velocities = [
    [velocity_wrapper(m, radii[i], θ[j], a, M, lum, correlation_length, mach) for j in 1:length(θ)]
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

keplerian_mags = [
    sqrt(sum(x -> x^2, keplerian_velocities[i][2:4]))
    for i in 1:length(radii)
]

difference_fracs = [
    [difference_mags[i][j] / keplerian_mags[i] for j in 1:length(θ)]
    for i in 1:length(radii)
]

# Find the minimum and maximum finite values in the nested array
min_val = minimum(filter(isfinite, collect(Iterators.flatten(difference_fracs))))
max_val = maximum(filter(isfinite, collect(Iterators.flatten(difference_fracs))))

for i in 1:length(radii)
    for j in 1:length(θ)
        if difference_fracs[i][j] == -Inf
            difference_fracs[i][j] = 2*min_val
        elseif difference_fracs[i][j] == Inf
            difference_fracs[i][j] = 2*max_val
        end
    end
end

difference_field = [100*difference_fracs[i][j] for i in 1:length(radii), j in 1:length(θ)]

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

# Define a custom color gradient with a logarithmic scale


begin
    plt = heatmap(
        θ,
        radii,
        difference_field,
        projection = :polar,
        title = "Magnitude of Velocity Difference", 
        colorbar_title = "Change in velocity compared to Keplerian (%)",
        grid = true,
        yticks = false,
        
    )
end