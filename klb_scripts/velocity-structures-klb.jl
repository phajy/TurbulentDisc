# Script includes the various velocity functions used to model turbulence

using Gradus, Plots, CoherentNoise, Statistics, StaticArrays, LinearAlgebra

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
function turb_random(m, r, correlation_length)
    keplerian = Gradus.CircularOrbits.fourvelocity(m, r)
    # add a random number between 0 and the sound speed, rescaled a little bit
    vt = SVector((1e-1 * randn() * soundspeed(r) for i in 1:4)...)

    v = keplerian .+ vt

    # now we need to ensure that the velocity has magnitude -1. This will
    # depend on the position `x`
    #x = SVector(0.0, r, π/2, 0.0)
    #Gradus.constrain_all(m, x, v, 0.0)
    v
end

# Perlin noise with finer variations and jitter
function turb_perlin(m, r, theta, correlation_length)
    
    keplerian = Gradus.CircularOrbits.fourvelocity(m, r)

    # Convert to cartesian coordinates
    x = r * cos(theta)
    y = r * sin(theta)

    scale_factor = 0.2

    # Generate Perlin noise
    perlin_noise = perlin_2d(seed=1)

    noise = scale_factor * sample(perlin_noise, x / correlation_length, y / correlation_length)

    # Add a small random jitter to break large-scale coherence/ smoothness
    noise += 0.01 * randn()

    vt = SVector(0, noise, 0, 0)
    
    # Add noise to Keplerian velocity
    v = keplerian .+ vt

    x = SVector(0.0, r, π/2, 0.0)

    # Ensure magnitude of 1
    Gradus.constrain_all(m, x, v, 1.0)
    
    return v
end

# fractional Brownian motion (fBm)
function turb_fbm(m, r, theta, correlation_length)
    
    keplerian = Gradus.CircularOrbits.fourvelocity(m, r)

    # fBm takes cartesian coordinates
    x = r * cos(theta)
    y = r * sin(theta)

    scale_factor = 0.2 

    # generate fBm noise 
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


function turbulent_structure(m, r, θ; type, correlation_length)
    if type == :random
        return turb_random(m, r, correlation_length)
    elseif type == :perlin
        return turb_perlin(m, r, θ, correlation_length)
    elseif type == :fbm
        return turb_fbm(m, r, θ, correlation_length)
    else
        throw(ArgumentError("Unknown turbulence type: $type"))
    end
end

# ---------- visualising the noise distributions (perlin against fBm), and how changing the scale factor affects these: ----------

# test params
m = KerrMetric(1.0, 0.998)  # Metric parameter (adjust as needed)
r_values = range(2.0, 10.0, length=100)  # Radial range
theta_values = range(0, 2π, length=100)  # Angular range
correlation_length = 2.0

# extrema vals
perlin_values = Float64[]
fbm_values = Float64[]

for r in r_values
    for θ in theta_values
        # Extract noise component (2nd index of vt)
        push!(perlin_values, turb_perlin(m, r, θ, correlation_length)[2])
        push!(fbm_values, turb_fbm(m, r, θ, correlation_length)[2])
    end
end

perlin_min, perlin_max = extrema(perlin_values)
fbm_min, fbm_max = extrema(fbm_values)

println("Perlin Noise Extrema: min = $perlin_min, max = $perlin_max")
println("fBm Noise Extrema: min = $fbm_min, max = $fbm_max")

histogram(perlin_values, bins=50, alpha=0.6, label="Perlin Noise", normalize=:pdf)
histogram!(fbm_values, bins=50, alpha=0.6, label="fBm Noise", normalize=:pdf, 
    title="Noise Distribution Comparison", xlabel="Velocity Perturbation", ylabel="Probability Density")



# ---------- visualising spatial structure of each noise functions (heatmap): ----------
perlin_matrix = reshape(perlin_values, (100, 100))
fbm_matrix = reshape(fbm_values, (100, 100))

heatmap(perlin_matrix, title="Perlin Noise Spatial Structure", xlabel="x", ylabel="y")
heatmap(fbm_matrix, title="fBm Noise Spatial Structure", xlabel="x", ylabel="y")
