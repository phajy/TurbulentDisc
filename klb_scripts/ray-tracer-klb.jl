# import the velocity functions
include("velocity-structures-klb.jl")

# Define correlation length
correlation_length = 1

function turbulent_redshift(metric, x_obs, vel_func, correlation_length)
    # Metric matrix at the observer's position
    g_obs = Gradus.metric(metric, x_obs)
    # Fixed stationary observer velocity
    v_obs = SVector{4, eltype(x_obs)}(1, 0, 0, 0)

    # Internal closure function
    function _internal_turbulent_redshift(m::AbstractMetric, gp, t)
        # Pass correlation_length to vel_func
        v_disc = vel_func(m, gp.x[2], gp.x[4], correlation_length)

        g = Gradus.metric(m, gp.x)
        Gradus.RedshiftFunctions._redshift_dotproduct(g, v_disc, g_obs, v_obs, gp)
    end

    return PointFunction(_internal_turbulent_redshift)
end

# Wrapper function that passes turb_func
function velocity_wrapper(m, r, theta, correlation_length)
    return turb_fbm(m, r, theta, correlation_length)
end

# Parameters and initial setup
m = KerrMetric(1.0, 0.998)
x = SVector(0.0, 1_000.0, deg2rad(80), 0.0)
d = ThinDisc(Gradus.isco(m), 50.0)

# Compose our turbulent redshift function with a filter function to remove points outside of the ISCO
redshift_pf = turbulent_redshift(m, x, velocity_wrapper, correlation_length)
pf = redshift_pf ∘ ConstPointFunctions.filter_intersected()

# Render geodesics using the modified point function 
α, β, img = rendergeodesics(
    m,
    x,
    d,
    20_000.0,  # Maximum integration time
    αlims = (-30, 30), 
    βlims = (-20, 20),
    image_width = 800,
    image_height = 400,
    verbose = true,
    pf = pf,
)

heatmap(α, β, img, aspect_ratio = 1)

# Create and plot line profile 
bins = collect(range(0.0, 2.0, 200))
plane = PolarPlane(GeometricGrid(); Nr = 1000, Nθ = 1000, r_max = 5 * d.outer_radius)
_, f = lineprofile(m, x, d, redshift_pf = pf, verbose = true, method = BinningMethod(), bins = bins, plane = plane)

# Set whatever is in the last bin to 0 as it's most likely a noise contribution
f[end] = 0
plot(bins, f, legend = false)

# next step: to force correlation_length to increase with r, simulating much a much more chaotic, granular turbulence close 
# to the black hole, and larger scale, less turbulent environment at larger r.