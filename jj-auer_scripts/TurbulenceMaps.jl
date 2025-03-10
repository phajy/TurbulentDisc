# -------- ONLY DEFINE FUNCTIONS AND STRUCTS IN THIS FILE ------------- #
# Use a seperate runner script (see e.g. plot-velocity.jl) to create
# an executable script. This file will be `include`d in other files, and
# any top-level code will be executed in all of them.

using Gradus, Plots, CoherentNoise
include("SoundandRadialSpeed.jl")

function logrange(first, last, num)
    10 .^ collect(range(log10(first), log10(last), num))
end

# -------- HERE IS WHERE YOU NEED TO NOW MODIFY THE VELOCITIES -------- #


function turbulence_random(m, r, correlation_length)
    keplerian = Gradus.CircularOrbits.fourvelocity(m, r)
    # add a random number between 0 and the sound speed, rescaled a little bit
    vt = SVector((1e-1 * randn() * c(m, r, 0.998, 1, 1e46) for i in 1:4)...)

    v = keplerian .+ vt

    # now we need to ensure that the velocity has magnitude -1. This will
    # depend on the position `x`
    # x = SVector(0.0, r, π/2, 0.0)
    # Gradus.constrain_all(m, x, v, 0.0)

    # NOTE: I have temporarily removed the normalisation, because for this
    # ad-hoc turbulent noise, it fails to constrain the velocity vector **and**
    # make the ray-traced image look good. I prioritise aesthetics over accuracy
    # here.
    v
end

function turbulence_perlin(m, r, theta, a, M, lum, correlation_length, mach; noise_r = perlin_2d(seed=4), noise_p = perlin_2d(seed=5), noise_a = perlin_2d(seed=6))

    f = inv(correlation_length)

    keplerian = Gradus.CircularOrbits.fourvelocity(m, r)

    # Convert polar to Cartesian coordinates
    x = r * cos(theta)
    y = r * sin(theta)
    
    # Generate Perlin noise based on the Cartesian coordinates
    noise_val_r = (1/sqrt(3))*mach*c(m, r, a, M, lum)*sample(noise_r, f*x, f*y)
    noise_val_p = (1/sqrt(3))*mach*c(m, r, a, M, lum)*sample(noise_p, f*x, f*y)
    noise_val_a = (1/sqrt(3))*mach*c(m, r, a, M, lum)*sample(noise_a, f*x, f*y)

    vt = SVector(0, noise_val_r, noise_val_p, noise_val_a)
    v = keplerian + vt

    # Ensure velocity constraint
    x = SVector(0.0, r, π/2, 0.0)
    Gradus.constrain_all(m, x, v, 1.0)

    return v
            
end

function turbulence_fbm(m, r, theta, a, M, ratio, alpha, correlation_length, mach; noise_r = fbm_fractal_2d(seed=1), noise_p = fbm_fractal_2d(seed=2), noise_a = fbm_fractal_2d(seed=3))

    f = inv(correlation_length)

    keplerian = Gradus.CircularOrbits.fourvelocity(m, r)

    # Convert polar to Cartesian coordinates
    x = r * cos(theta)
    y = r * sin(theta)
    
    # Generate Perlin noise based on the Cartesian coordinates
    noise_val_r = (1/sqrt(3))*mach*SoundSpeed(m, r, a, M, ratio)*sample(noise_r, f*x, f*y)
    noise_val_p = (1/sqrt(3))*mach*SoundSpeed(m, r, a, M, ratio)*sample(noise_p, f*x, f*y)
    noise_val_a = (1/sqrt(3))*mach*SoundSpeed(m, r, a, M, ratio)*sample(noise_a, f*x, f*y)

    radial_inflow = RadialSpeed(m, r, a, alpha, M, ratio)

    v_additional = SVector(0, noise_val_r + radial_inflow, noise_val_p, noise_val_a)
    v = keplerian + v_additional

    # Ensure velocity constraint
    x = SVector(0.0, r, π/2, 0.0)
    Gradus.constrain_all(m, x, v, 1.0)

    return v
            
end
