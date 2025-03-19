"""
This script defines velocity functions for turbulent accretion disc models using Perlin and fractional Brownian motion (fBm) noise.
 It computes the effective velocity of fluid elements in a relativistic accretion disc, incorporating:
 - A Keplerian velocity component from the Gradus CircularOrbits module.
 - A turbulent velocity perturbation based on either Perlin noise or fBm noise.
 - A radial inflow component derived from an external radial velocity function.
"""

using Gradus, Plots, CoherentNoise
include("updated-sound-speed.jl")

function turbulence_perlin(m, r, theta, a, M, ratio, alpha, correlation_length, mach; noise_r = perlin_2d(seed=1), noise_p = perlin_2d(seed=2), noise_a = perlin_2d(seed=3))

    v_keplerian = Gradus.CircularOrbits.fourvelocity(m, r)

    # Convert polar to Cartesian coordinates
    x = r * cos(theta)
    y = r * sin(theta)
    
    # Generate Perlin noise based on the Cartesian coordinates
    noise_val_r = (1/sqrt(3)) * mach * SoundSpeed(m, r, a, M, ratio) * sample(noise_r, x/correlation_length, y/correlation_length)
    noise_val_p = (1/sqrt(3)) * mach * SoundSpeed(m, r, a, M, ratio) * sample(noise_p, x/correlation_length, y/correlation_length)
    noise_val_a = (1/sqrt(3)) * mach * SoundSpeed(m, r, a, M, ratio) * sample(noise_a, x/correlation_length, y/correlation_length)

    radial_inflow = RadialSpeed(m, r, a, alpha, M, ratio)

    v_turb = SVector(0, noise_val_r + radial_inflow, noise_val_p, noise_val_a)
    v_eff = v_keplerian + v_turb

    # Ensure velocity constraint
    x = SVector(0.0, r, π/2, 0.0)
    Gradus.constrain_all(m, x, v_eff, 1.0)

    return v_eff
            
end

function turbulence_fbm(m, r, theta, a, M, ratio, alpha, correlation_length, mach; noise_r = fbm_fractal_2d(seed=4), noise_p = fbm_fractal_2d(seed=5), noise_a = fbm_fractal_2d(seed=6))

    v_keplerian = Gradus.CircularOrbits.fourvelocity(m, r)

    # Convert polar to Cartesian coordinates
    x = r * cos(theta)
    y = r * sin(theta)
    
    # Generate Perlin noise based on the Cartesian coordinates
    noise_val_r = (1/sqrt(3)) * mach * SoundSpeed(m, r, a, M, ratio) * sample(noise_r, x/correlation_length, y/correlation_length)
    noise_val_p = (1/sqrt(3)) * mach * SoundSpeed(m, r, a, M, ratio) * sample(noise_p, x/correlation_length, y/correlation_length)
    noise_val_a = (1/sqrt(3)) * mach * SoundSpeed(m, r, a, M, ratio) * sample(noise_a, x/correlation_length, y/correlation_length)

    radial_inflow = RadialSpeed(m, r, a, alpha, M, ratio)

    v_turb = SVector(0, noise_val_r + radial_inflow, noise_val_p, noise_val_a)
    v_eff = v_keplerian + v_turb

    # Ensure velocity constraint
    x = SVector(0.0, r, π/2, 0.0)
    Gradus.constrain_all(m, x, v_eff, 1.0)

    return v_eff
            
end