# import functions
include("VelocityStructures.jl")

m = KerrMetric(1.0, 0.998)

# log range from the innermost stable circular orbit (inner edge of the
# accretion disc) out to some arbitrary radius
radii = logrange(Gradus.isco(m), 1000, 200)

# if you need it, these will be the azimuthal angles to each "radial" bin
θ = collect(range(0, 2π, 200))

keplerian_velocities = Gradus.CircularOrbits.fourvelocity.(m, radii)

# then apply however your velocity function works to the keplerian_velocities
turbulent_velocities = turbulence_random.(m, radii, θ)

# --------------------------------------------------------------------- #

# some visualisation methods (use the ones that make the most sense to you)
# note these are all only using the 4th component of the velocity vector (which
# is the azimuthal component):
# 1 - time
# 2 - radial (r)
# 3 - poloidal (θ)
# 4 - azimuthal (ϕ)
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
    turbulent_field = [log10(turbulent_structure(m, r, 0.0)[4]) for r in radii, angle in θ]
    heatmap(θ, radii, turbulent_field, projection = :polar, title = "turbulent")
end

begin
    comparison_map = [
        (turbulent_structure(m, r, 0.0)[4] - v[4]) / v[4] for
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
