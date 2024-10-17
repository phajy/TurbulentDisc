using Gradus, Plots, Noise

rng = MersenneTwister(1)

function logrange(first, last, num)
    10 .^ collect(range(log10(first), log10(last), num))
end

"""
This function is an order of magnitude approximation for the Novikov-Thorne
temperature profile for a radiation-pressure dominated accretion disc:

    T ∝ R^(-3/4)

And the ideal gas temperature and sound speed relation

    cₛ^2 ∝ T^4

!!! note

    There will be some dependence on the metric in this function, but to get
    going this is probably good enough.
"""
function soundspeed(r)
    T = r^(-3 / 4)
    T^2
end

m = KerrMetric(1.0, 0.998)

# log range from the innermost stable circular orbit (inner edge of the
# accretion disc) out to some arbitrary radius
radii = logrange(Gradus.isco(m), 1000, 200)

# if you need it, these will be the azimuthal angles to each "radial" bin
thetas = collect(range(0, 2π, 100))

keplerian_velocities = Gradus.CircularOrbits.fourvelocity.(m, radii)

# -------- HERE IS WHERE YOU NEED TO NOW MODIFY THE VELOCITIES -------- #

# below is just a trivial example
function turbulent_test(m, r, correlation_length)

    Plots.Random.seed!(1)

    keplerian = Gradus.CircularOrbits.fourvelocity(m, r)
    # add a random number between 0 and the sound speed, rescaled a little bit
    vt = SVector((2e-1 * randn() * soundspeed(r) for _ = 1:4)...)

    v = keplerian .+ vt

    # now we need to ensure that the velocity has magnitude -1. This will
    # depend on the position `x`
    x = SVector(0.0, r, 0.0, 0.0)
    Gradus.constrain_all(m, x, v, 1.0)
end

function turbulent_perlin(m, radii, thetas)

    keplerian_map = [Gradus.CircularOrbits.fourvelocity(m, r) for r in radii]
    turbulence_map = Matrix{SVector{4, Float64}}(undef, length(radii), length(thetas))

    # Generate Perlin noise for each (r, θ) pair
    perlin = perlin_2d()

    for i in 1:length(radii)
        for j in 1:length(thetas)
            # Convert polar to Cartesian coordinates
            r = radii[i]
            θ = thetas[j]
            x = r * cos(θ)
            y = r * sin(θ)
            
            # Generate Perlin noise based on the Cartesian coordinates
            noise_val = noise(perlin_noise, x, y)

            v = keplerian_map[i] .+ SVector(0, 0, 0, noise_val)

            x = SVector(0.0, r, 0.0, 0.0)
            turbulence_map[i, j] = Gradus.constrain_all(m, x, v, 1.0)

        end
    end

    return turbulence_map
            
end

# then apply however your velocity function works to the keplerian_velocities
#turbulent_velocities = turbulent_perlin(m, radii, thetas)

turbulent_velocities = turbulent_perlin(m, radii, thetas)


# --------------------------------------------------------------------- #

# some visualisation methods (use the ones that make the most sense to you)
# note these are all only using the 4th component of the velocity vector (which
# is the azimuthal component):
# 1 - time
# 2 - radial (r)
# 3 - poloidal (thetas)
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
        title = "(turbulent - keplerian) / keplerian",
    )
end