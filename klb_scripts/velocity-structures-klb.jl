# Script includes the various velocity functions used to model turbulence

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

# trivial, random example:
function turb_rand(m, r, correlation_length)
    keplerian = Gradus.CircularOrbits.fourvelocity(m, r)
    # selects random number and scaled by sound speed as a function of r
    vt = SVector((1e-1 * randn() * soundspeed(r) for i in 1:4)...)

    # add random perturbations to Keplerian velocities
    v = keplerian .+ vt

    # now we need to ensure that the velocity has magnitude -1. This will
    # depend on the position `x`
    x = SVector(0.0, r, π/2, 0.0)

    Gradus.constrain_all(m, x, v, 1.0)
end

# Perlin noise
function turb_perlin(m, r, theta, correlation_length)
    
    keplerian = Gradus.CircularOrbits.fourvelocity(m, r)

    # perlin_noise takes cartesian coordinates
    x = r * cos(theta)
    y = r * sin(theta)

    scale_factor = 0.2 

    # generate Perlin noise for x,y sample
    perlin_noise = perlin_2d(seed=1)

    # sample noise at scaled coordinates (x, y)
    noise = scale_factor * sample(perlin_noise, x / correlation_length, y / correlation_length)

    vt = SVector(0, noise, 0, 0)

    # add noise to Keplerian velocities
    v = keplerian .+ vt

    x = SVector(0.0, r, π/2, 0.0)

    # ensure magnitude of 1
    Gradus.constrain_all(m, x, v, 1.0)

end

# fractional Brownian motion (fBm)
function turb_fbm(m, r, theta, correlation_length)
    
    keplerian = Gradus.CircularOrbits.fourvelocity(m, r)

    # fBm takes cartesian coordinates
    x = r * cos(theta)
    y = r * sin(theta)

    scale_factor = 0.2 

    # generate fBm noise with a specified parameters
    fbm_noise = fbm_fractal_2d(seed=1, octaves=4, frequency=1.0, lacunarity=2.0, persistence=0.5)

    # sample the fBm noise at scaled coordinates (x, y) 
    noise = scale_factor * sample(fbm_noise, x / correlation_length, y / correlation_length)

    vt = SVector(0, noise, 0, 0)

    # add the noise to the Keplerian velocity
    v = keplerian .+ vt

    x = SVector(0.0, r, π/2, 0.0)

    # constrain the magnitude to 1
    Gradus.constrain_all(m, x, v, 1.0)

    return v
end