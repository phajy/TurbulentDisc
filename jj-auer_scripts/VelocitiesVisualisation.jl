# import functions
include("TurbulenceMaps.jl")

m = KerrMetric(1.0, 0.998)

# log range from the innermost stable circular orbit (inner edge of the
# accretion disc) out to some arbitrary radius
radii = logrange(Gradus.isco(m), 1000, 200)

# if you need it, these will be the azimuthal angles to each "radial" bin
θ = collect(range(0, 2π, 200))

keplerian_velocities = Gradus.CircularOrbits.fourvelocity.(m, radii)

a=0.998
M=1.0
lum=7.2e42
correlation_length=1
mach=1

# then apply however your velocity function works to the keplerian_velocities
turbulent_velocities = turbulence_fbmfractal(m, radii, θ, a, M, lum, correlation_length, mach)

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
        [i[4] for i in keplerian_velocities],
        xlabel = "r",
        ylabel = "vᵩ",
        xscale = :log10,
        yscale = :log10,
        label = "keplerian",
    )
    plot!(radii, [i[4] for i in turbulent_velocities], label = "turbulent")
end


begin
    keplerian_field =
        [log10(i[4]) for (i, r) in zip(keplerian_velocities, radii), angle in θ]
    heatmap(θ, radii, keplerian_field, projection = :polar, title = "keplerian")
end

begin
    turbulent_field = [log10(turbulence_perlin(m, r, 0.0)[4]) for r in radii, angle in θ]
    heatmap(θ, radii, turbulent_field, projection = :polar, title = "turbulent")
end
"""


begin
    comparison_map = [
        (turbulence_fbmfractal(m, r, 0.0, 0.998, 1.0, 7.2e42, 1, 1)[4] - v[4]) / v[4] for
        (v, r) in zip(keplerian_velocities, radii), angle in θ
    ]
    heatmap(
        θ,
        radii,
        comparison_map,
        projection = :polar,
        title = "(Turbulent - Keplerian) / Keplerian",
    )
end
