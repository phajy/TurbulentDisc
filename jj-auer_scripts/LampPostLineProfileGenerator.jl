using Gradus

include("TurbulenceMaps.jl")

# to include our turbulent velocity map, we need to specify a new
# point function that calculates the redshift with our custom velocity

m = KerrMetric(1.0, 0.998)
x = SVector(0.0, 1_000.0, deg2rad(80), 0.0)
d = ThinDisc(Gradus.isco(m), 50.0)

function turbulent_redshift(metric, x_obs, vel_func)
    # metric matrix at the observer's position
    g_obs = Gradus.metric(metric, x_obs)
    # fixed stationary observer velocity
    v_obs = SVector{4,eltype(x)}(1, 0, 0, 0)

    # internal closure
    function _internal_turbulent_redshift(m::AbstractMetric, gp, t)
        v_disc = turbulence_fbmfractal(m, gp.x[2], gp.x[4])

        g = Gradus.metric(m, gp.x)
        Gradus.RedshiftFunctions._redshift_dotproduct(g, v_disc, g_obs, v_obs, gp)
    end

    return PointFunction(_internal_turbulent_redshift)
end


# NOTE: this function is not valid within the ISCO, so we need to make sure that
# we **always** set the inner radius of the disc to the ISCO
function velocity_wrapper(m, r, phis)
    return turbulence_fbmfractal(m, r, phis)
end

# compose our turbulent redshift function with a filter function to remove those
# points outside of the ISCO
redshift_pf = turbulent_redshift(m, x, velocity_wrapper)
pf = redshift_pf ∘ ConstPointFunctions.filter_intersected()

model = LampPostModel(h = 10.0)
ε(r) = emissivity_profile(m, d, model)

bins = collect(range(0.1, 1.4, 200))

function calculate_line_profile(m, x, d, bins)
    ε(r) = r^(-4)  # Define emissivity function with given index q
    _, f = lineprofile(
        m,
        x,
        d,
        method = BinningMethod(),
        callback = domain_upper_hemisphere(),
        verbose = true,
        bins = bins
    )
    return f
end

flux = calculate_line_profile(m, x, d, bins)

# Plot the line profile with annotations for inclination angle and emissivity index
plot(
    bins, flux,
    xlabel = "Redshift",
    ylabel = "Flux (arbitrary units)",
    legend = false,
)
