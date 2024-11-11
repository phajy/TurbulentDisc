# -------- ONLY DEFINE FUNCTIONS AND STRUCTS IN THIS FILE ------------- #
# Use a seperate runner script (see e.g. plot-velocity.jl) to create
# an executable script. This file will be `include`d in other files, and
# any top-level code will be executed in all of them.

using Gradus, Plots, CoherentNoise

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

# -------- HERE IS WHERE YOU NEED TO NOW MODIFY THE VELOCITIES -------- #


function turbulence_random(m, r, correlation_length)
    keplerian = Gradus.CircularOrbits.fourvelocity(m, r)
    # add a random number between 0 and the sound speed, rescaled a little bit
    vt = SVector((1e-1 * randn() * soundspeed(r) for i in 1:4)...)

    v = keplerian .+ vt

    # now we need to ensure that the velocity has magnitude -1. This will
    # depend on the position `x`
    # x = SVector(0.0, r, π/2, 0.0)
    # Gradus.constrain_all(m, x, v, 0.0)

    # NOTE: I have temporarily removed the normalisation, because for this
    # ad-hoc turbulent noise, it fails to constrain the velocity vector **and**
    # make the ray-traced image look good. I prioritise aesthetics over accuracy
    # here.
    Gradus.constrain_all(m, x, v, 1.0)
end


perlin_noise = perlin_2d(seed=1)

function turbulence_perlin(m, r, theta)

    intensity=0.2

    keplerian = Gradus.CircularOrbits.fourvelocity(m, r)

    # Convert polar to Cartesian coordinates
    x = r * cos(theta)
    y = r * sin(theta)

    
    # Generate Perlin noise based on the Cartesian coordinates
    noise_val = intensity*sample(perlin_noise, x, y)

    # @show noise_val

    vt = SVector(0, noise_val, 0, 0)
    v = keplerian + vt
    x_disc = SVector(0.0, r, π/2, 0.0)
    Gradus.constrain_all(m, x_disc, v, 1.0)
end
