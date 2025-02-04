# Script includes the various velocity functions used to model turbulence

using Gradus, Plots, CoherentNoise, Statistics, StaticArrays, LinearAlgebra
include("pariev-bromley-equations-klb.jl")  # Import sound speed function

# Wrapper function for sound speed using Pariev & Bromley prescription
c(m, r, a, M, L, L_edd, r_ms, epsilon) = sound_speed_ratio(r, a, epsilon, L, L_edd, r_ms, M)

function logrange(first, last, num)
    10 .^ collect(range(log10(first), log10(last), num))
end

# -------- TURBULENCE FUNCTIONS WITH SOUND SPEED NORMALIZATION -------- #

# Random turbulence
function turb_random(m, r, a, M, L, L_edd, r_ms, epsilon, correlation_length)
    keplerian = Gradus.CircularOrbits.fourvelocity(m, r)

    # Generate random perturbation normalized to sound speed
    vt = SVector((1e-1 * randn() * c(m, r, a, M, L, L_edd, r_ms, epsilon) for i in 1:4)...)

    v = keplerian .+ vt
    return v
end

# Perlin noise turbulence
function turb_perlin(m, r, theta, a, M, L, L_edd, r_ms, epsilon, correlation_length, mach)
    keplerian = Gradus.CircularOrbits.fourvelocity(m, r)

    # Convert to Cartesian coordinates
    x = r * cos(theta)
    y = r * sin(theta)

    scale_factor = 1

    # Generate Perlin noise
    perlin_noise = perlin_2d(seed=1)

    # Sample Perlin noise 
    noise = scale_factor * sample(perlin_noise, x / correlation_length, y / correlation_length)
    noise *= mach * c(m, r, a, M, L, L_edd, r_ms, epsilon) # normalise to sound speed

    vt = SVector(0, noise, 0, 0)
    
    # Add noise to Keplerian velocity
    v = keplerian .+ vt

    # Ensure velocity constraint
    x = SVector(0.0, r, π/2, 0.0)
    Gradus.constrain_all(m, x, v, 1.0)

    return v
end

# Fractional Brownian motion (fBm) turbulence
function turb_fbm(m, r, theta, a, M, L, L_edd, r_ms, epsilon, correlation_length, mach)
    keplerian = Gradus.CircularOrbits.fourvelocity(m, r)

    # Convert to Cartesian coordinates
    x = r * cos(theta)
    y = r * sin(theta)

    scale_factor = 1

    # Generate fBm noise
    fbm_noise = fbm_fractal_2d(seed=1, octaves=4, frequency=1.0, lacunarity=2.0, persistence=0.5)

    # Sample fBm noise 
    noise = scale_factor * sample(fbm_noise, x / correlation_length, y / correlation_length)
    noise *= mach * c(m, r, a, M, L, L_edd, r_ms, epsilon) # normalise to sound speed

    vt = SVector(0, noise, 0, 0)

    # Add noise to Keplerian velocity
    v = keplerian .+ vt

    # Ensure velocity constraint
    x = SVector(0.0, r, π/2, 0.0)
    Gradus.constrain_all(m, x, v, 1.0)

    return v
end

function turbulent_structure(m, r, θ; type, a, M, L, L_edd, r_ms, epsilon, correlation_length, mach)
    if type == :random
        return turb_random(m, r, a, M, L, L_edd, r_ms, epsilon, correlation_length)
    elseif type == :perlin
        return turb_perlin(m, r, θ, a, M, L, L_edd, r_ms, epsilon, correlation_length, mach)
    elseif type == :fbm
        return turb_fbm(m, r, θ, a, M, L, L_edd, r_ms, epsilon, correlation_length, mach)
    else
        throw(ArgumentError("Unknown turbulence type: $type"))
    end
end

# ---------- VISUALISATION CODE ---------- #

# Test parameters
m = KerrMetric(1.0, 0.998)  # Metric parameter
r_values = range(2.0, 10.0, length=100)  # Radial range
theta_values = range(0, 2π, length=100)  # Angular range
correlation_length = 2.0
a = 0.998  # Black hole spin
M = 1.0  # Black hole mass in geometrized units
L = 1e46  # Luminosity in ergs/s
L_edd = eddington_luminosity(M)  # Compute Eddington luminosity
r_ms = Gradus.isco(m)  # Compute ISCO
epsilon = 0.1  # Efficiency factor
mach = 1 # Mach number

# Collect values
perlin_values = Float64[]
fbm_values = Float64[]

for r in r_values
    for θ in theta_values
        push!(perlin_values, turb_perlin(m, r, θ, a, M, L, L_edd, r_ms, epsilon, correlation_length, mach)[2])
        push!(fbm_values, turb_fbm(m, r, θ, a, M, L, L_edd, r_ms, epsilon, correlation_length, mach)[2])
    end
end

perlin_min, perlin_max = extrema(perlin_values)
fbm_min, fbm_max = extrema(fbm_values)

println("Perlin Noise Extrema: min = $perlin_min, max = $perlin_max")
println("fBm Noise Extrema: min = $fbm_min, max = $fbm_max")

histogram(perlin_values, bins=50, alpha=0.6, label="Perlin Noise", normalize=:pdf)
histogram!(fbm_values, bins=50, alpha=0.6, label="fBm Noise", normalize=:pdf, 
    title="Noise Distribution Comparison", xlabel="Velocity Perturbation", ylabel="Probability Density")

# ---------- HEATMAP VISUALIZATION ---------- #

perlin_matrix = reshape(perlin_values, (100, 100))
fbm_matrix = reshape(fbm_values, (100, 100))

heatmap(perlin_matrix, title="Perlin Noise Spatial Structure", xlabel="x", ylabel="y")
heatmap(fbm_matrix, title="fBm Noise Spatial Structure", xlabel="x", ylabel="y")
