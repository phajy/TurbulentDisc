include("VelocityStructures.jl")
include("velocity-structures.jl")

# to include our turbulent velocity map, we need to specify a new
# point function that calculates the redshift with our custom velocity

"""
    turbulent_redshift(metric, x_obs, vel_func)

A [`PointFunction`](@ref) internal that uses the custom velocity map `vel_func`
to create a redshift point function.

`vel_func` should be a function with the following signature

```julia
function turbulent_structure(m, r::Number, ϕ::Number)::SVector{4}
```
"""
function turbulent_redshift(metric, x_obs, vel_func)
    # metric matrix at the observer's position
    g_obs = Gradus.metric(metric, x_obs)
    # fixed stationary observer velocity
    v_obs = SVector{4,eltype(x)}(1, 0, 0, 0)

    # internal closure
    function _internal_turbulent_redshift(m::AbstractMetric, gp, t)
        v_disc = vel_func(m, gp.x[2], gp.x[4])

        g = Gradus.metric(m, gp.x)
        Gradus.RedshiftFunctions._redshift_dotproduct(g, v_disc, g_obs, v_obs, gp)
    end

    return PointFunction(_internal_turbulent_redshift)
end

# NOTE: this function is not valid within the ISCO, so we need to make sure that
# we **always** set the inner radius of the disc to the ISCO
function velocity_wrapper(m, r, phis)
    return turbulence_perlin(m, r, phis)
end



m = KerrMetric(1.0, 0.998)
x = SVector(0.0, 1_000.0, deg2rad(80), 0.0)
d = ThinDisc(Gradus.isco(m), 50.0)

# compose our turbulent redshift function with a filter function to remove those
# points outside of the ISCO
redshift_pf = turbulent_redshift(m, x, velocity_wrapper)
pf = redshift_pf ∘ ConstPointFunctions.filter_intersected()



α, β, img = rendergeodesics(
    m,
    x,
    d,
    # maximum integration time
    20_000.0,
    αlims = (-30, 30), 
    βlims = (-20, 20),
    image_width = 800,
    image_height = 400,
    verbose = true,
    pf = pf,
)

heatmap(α, β, img, aspect_ratio = 1)

#=

# to create a lineprofile, probably best to use the binning method until some
# correlation length exists, else the lineprofiles will amplify noise to due
# discontinuous redshift functions
bins = collect(range(0.0, 2.0, 200))
_, f = lineprofile(m, x, d, redshift_pf = pf, verbose = true, method = BinningMethod(), bins = bins)

# set whatever is in the last bin to 0 as it's most likely a noise contribution
f[end] = 0
plot(bins, f)
=#