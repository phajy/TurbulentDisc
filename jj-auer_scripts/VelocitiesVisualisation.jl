# import functions
using Gradus, Plots
include("TurbulenceMaps.jl")

bins=200

a = 0.80
M = 1.0
lum = 7.2e42
correlation_length = 5
mach = 20000
outer_radius = 40.0
m = KerrMetric(M, a)
d = Gradus.ThinDisc(Gradus.isco(m), outer_radius)
incl_angle = 80
x = SVector(0.0, 1_000.0, deg2rad(incl_angle), 0.0)

radii = logrange(Gradus.isco(m), outer_radius, bins)
θ = collect(range(0, 2π, bins))

function velocity_wrapper(m, r, phis, a, M, lum, correlation_length, mach)
    return turbulence_fbm(m, r, phis, a, M, lum, correlation_length, mach)
end

function turbulent_redshift(metric, x_obs, vel_func, a, M, lum, correlation_length, mach)
    # metric matrix at the observer's position
    g_obs = Gradus.metric(metric, x_obs)
    # fixed stationary observer velocity
    v_obs = SVector{4,eltype(x)}(1, 0, 0, 0)

    # internal closure
    function _internal_turbulent_redshift(m::AbstractMetric, gp, t)
        v_disc = vel_func(m, gp.x[2], gp.x[4], a, M, lum, correlation_length, mach)

        g = Gradus.metric(m, gp.x)
        Gradus.RedshiftFunctions._redshift_dotproduct(g, v_disc, g_obs, v_obs, gp)
    end

    return PointFunction(_internal_turbulent_redshift)
end


# Collecting the velocity differences for the polar heatmap
begin

    keplerian_velocities = Gradus.CircularOrbits.fourvelocity.(m, radii)

    # Calculate turbulent velocities for each combination of r and θ
    turbulent_velocities = [
        [velocity_wrapper(m, radii[i], θ[j], a, M, lum, correlation_length, mach) for j in 1:length(θ)]
        for i in 1:length(radii)
    ]

    turbulent_velocities_slice = [turbulent_velocities[i][1] for i in 1:length(radii)]

    # Calculate the difference between turbulent and keplerian velocities
    velocity_differences = [
        [turbulent_velocities[i][j] - keplerian_velocities[i] for j in 1:length(θ)]
        for i in 1:length(radii)
    ]

    # Calculate magnitudes of the differences (excluding the first component)
    difference_mags = [
        [sqrt(sum(x -> x^2, velocity_differences[i][j][2:4])) for j in 1:length(θ)]
        for i in 1:length(radii)
    ]

    keplerian_mags = [
        sqrt(sum(x -> x^2, keplerian_velocities[i][2:4]))
        for i in 1:length(radii)
    ]

    difference_fracs = [
        [difference_mags[i][j] / keplerian_mags[i] for j in 1:length(θ)]
        for i in 1:length(radii)
    ]

    # Find the minimum and maximum finite values in the nested array
    min_val = minimum(filter(isfinite, collect(Iterators.flatten(difference_fracs))))
    max_val = maximum(filter(isfinite, collect(Iterators.flatten(difference_fracs))))

    for i in 1:length(radii)
        for j in 1:length(θ)
            if difference_fracs[i][j] == -Inf
                difference_fracs[i][j] = 2*min_val
            elseif difference_fracs[i][j] == Inf
                difference_fracs[i][j] = 2*max_val
            end
        end
    end

    difference_field = [100*difference_fracs[i][j] for i in 1:length(radii), j in 1:length(θ)]

end

# 1 - time
# 2 - radial (r)
# 3 - poloidal (θ)
# 4 - azimuthal (ϕ)

# MAKE SURE UNITS ARE CORRECT AND PLOTTED AT THE AXES - SHOULD BE RG


# 1D Plot of turbulence added

begin

    plt = plot(
    xlabel = "r (M))",
    ylabel = "Δv (unitless)",
    legend = :topright,
    )

    radial_differences = [(turbulent_velocities_slice[i][2] - keplerian_velocities[i][2]) for i in 1:bins]

    plot!(
        radii,
        radial_differences,
        label = "Radial Turbulence Added",
    )

    poloidal_differences = [(turbulent_velocities_slice[i][3] - keplerian_velocities[i][3]) for i in 1:bins]

    plot!(
        radii,
        poloidal_differences,
        label = "Poloidal Turbulence Added",
    )

    azimuthal_differences = [(turbulent_velocities_slice[i][4] - keplerian_velocities[i][4]) for i in 1:bins]

    plot!(
        radii,
        azimuthal_differences,
        label = "Azimuthal Turbulence Added",
    )

    display(plt)

end


# Turbulence polar heatmap

begin
    plt = heatmap(
        θ,
        #log10.(radii),
        radii,
        difference_field,
        projection = :polar,
        title = "Magnitude of Velocity Difference", 
        colorbar_title = "Change in velocity compared to Keplerian (%)",
        labelpad = 5,
        grid = true,
        yticks = false,

    )

    display(plt)

end


# Redshift heatmap

begin
    redshift_pf = turbulent_redshift(m, x, velocity_wrapper, a, M, lum, correlation_length, mach)
    print(redshift_pf)

    pf = redshift_pf ∘ ConstPointFunctions.filter_intersected()

    α, β, img = rendergeodesics(
        m,
        x,
        d,
        # maximum integration time
        2000.0,
        αlims = (-50, 50), 
        βlims = (-18, 20),
        image_width = 800,
        image_height = 400,
        verbose = true,
        pf = pf,
    )

    heatmap(α, β, img, aspect_ratio = 1, color=:inferno)
end



#savefig(plt, "Other/Figs/FillInHere.pdf")

