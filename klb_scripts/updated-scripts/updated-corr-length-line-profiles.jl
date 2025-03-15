using Plots, Gradus, Measures, LaTeXStrings
include("updated-velocity-functions.jl")
include("updated-sound-speed.jl")       

plot_font = "Computer Modern"
default(fontfamily=plot_font,
        linewidth=0.8, framestyle=:box, label=nothing, grid=false)
scalefontsizes(1)

M = 1.0                  # Black hole mass
a = 0.998                # Black hole spin parameter
ratio_thin = 0.3         # Eddington ratio for thin disc
ratio_thick = 0.5        # Eddington ratio for thick disc
h = 10.0                 # Lamppost corona height
m = KerrMetric(M, a)     # Kerr metric
inner_radius = Gradus.isco(m)  # Inner disc radius
outer_radius = 50.0      # Outer disc radius for thick disc
h_scale = 0.3            # Scale height for thick disc
bins = collect(range(0.1, 2.0, 200))
correlation_lengths = [0.5, 1, 10]  
inclination_angles = [30, 60, 75]
turbulence_model = "fbm" # 'perlin' or 'fbm'
mach = 1

opacity_values = Dict(0.5 => 1.0, 1 => 0.8, 10 => 0.6)
linewidth = 0.8

output_dir = "C:\\Users\\Kate\\project\\TurbulentDisc\\updated_klb_plots\\line_profiles\\turbulent\\correlation_length"
mkpath(output_dir)


# Function for turbulent redshift
function turbulent_redshift(m, x_obs, vel_func, a, M, ratio, alpha, correlation_length, mach)
    g_obs = Gradus.metric(m, x_obs)  
    v_obs = SVector{4, eltype(x_obs)}(1, 0, 0, 0) 

    function _internal_turbulent_redshift(m::AbstractMetric, gp, t)
        v_disc = vel_func(m, gp.x[2], gp.x[4], a, M, ratio, alpha, correlation_length, mach)
        g = Gradus.metric(m, gp.x)
        Gradus.RedshiftFunctions._redshift_dotproduct(g, v_disc, g_obs, v_obs, gp)
    end

    return PointFunction(_internal_turbulent_redshift)
end

# Wrapper for velocity function
function velocity_wrapper(m, r, phis, a, M, ratio, alpha, correlation_length, mach)
    if turbulence_model == "perlin"
        return turbulence_perlin(m, r, phis, a, M, ratio, alpha, correlation_length, mach)
    elseif turbulence_model == "fbm"
        return turbulence_fbm(m, r, phis, a, M, ratio, alpha, correlation_length, mach)
    end
end

# --- Function for Turbulent Line Profile ---
function calculate_turbulent_line_profile(m, x, d, bins, correlation_length, a, M, ratio, alpha, mach, ε)
    redshift_pf = turbulent_redshift(m, x, velocity_wrapper, a, M, ratio, alpha, correlation_length, mach)
    pf = redshift_pf ∘ ConstPointFunctions.filter_intersected()
    plane = PolarPlane(GeometricGrid(); Nr = 1000, Nθ = 1000, r_max = outer_radius)
    _, f = lineprofile(m, x, d, ε; bins, redshift_pf = pf, method = BinningMethod(), plane = plane)
    f[end] = 0
    return f
end

# --- Function for Zero-Turbulence (Laminar) Line Profile ---
function calculate_zero_turbulence_line_profile(m, x, d, bins, ε)
    _, f = lineprofile(m, x, d, ε; bins, method = BinningMethod(), callback = domain_upper_hemisphere(), verbose = true)
    f[end] = 0
    return f
end

function height_profile_exponential(ρ)
    if ρ <= inner_radius || ρ > outer_radius
        return 0  
    end
    return h_scale * ρ * exp(-ρ / outer_radius)
end

for incl_angle in inclination_angles
    x = SVector(0.0, 1_000.0, deg2rad(incl_angle), 0.0)

    # Thin Disc Plot
    plt_thin = plot(
        xlabel = L"\nu/\nu_e",
        ylabel = L"\textrm{Flux \ (Arbitrary \ Units)}",
        legend = :topleft,
        left_margin = [5mm 0mm], right_margin = [5mm 0mm],
        top_margin = [5mm 0mm], bottom_margin = [5mm 0mm],
        xlims = (0, 2.0)
    )

    m = KerrMetric(M, a)
    d_thin = ShakuraSunyaev(m, eddington_ratio=ratio_thin)
    ε_thin = emissivity_profile(m, d_thin, LampPostModel(h=h))

    f_thin_laminar = calculate_zero_turbulence_line_profile(m, x, d_thin, bins, ε_thin)
    plot!(plt_thin, bins, f_thin_laminar, label="Zero Turbulence", color=:black, alpha=1.0, linestyle=:solid, linewidth=linewidth)

    for correlation_length in correlation_lengths
        opacity = opacity_values[correlation_length] 
        f_thin_turb = calculate_turbulent_line_profile(m, x, d_thin, bins, correlation_length, a, M, ratio_thin, 0.1, mach, ε_thin)
        plot!(plt_thin, bins, f_thin_turb, label=L"l_{c} = %$correlation_length", color=:black, alpha=opacity, linestyle=:dash, linewidth=linewidth)
    end

    filename_thin = joinpath(output_dir, "$(turbulence_model)_line_profile_thin_incl$(incl_angle).pdf")
    savefig(plt_thin, filename_thin)
    display(plt_thin)

    # Thick Disc Plot
    plt_thick = plot(
        xlabel = L"\nu/\nu_e",
        ylabel = L"\textrm{Flux \ (Arbitrary \ Units)}",
        legend = :topleft,
        left_margin = [5mm 0mm], right_margin = [5mm 0mm],
        top_margin = [5mm 0mm], bottom_margin = [5mm 0mm],
        xlims = (0, 2.0)
    )

    d_thick = ThickDisc(height_profile_exponential)
    ε_thick = emissivity_profile(m, d_thick, LampPostModel(h=h))

    f_thick_laminar = calculate_zero_turbulence_line_profile(m, x, d_thick, bins, ε_thick)
    plot!(plt_thick, bins, f_thick_laminar, label="Zero Turbulence", color=:black, alpha=1.0, linestyle=:solid, linewidth=linewidth)

    for correlation_length in correlation_lengths
        opacity = opacity_values[correlation_length] 
        f_thick_turb = calculate_turbulent_line_profile(m, x, d_thick, bins, correlation_length, a, M, ratio_thick, 0.1, mach, ε_thick)
        plot!(plt_thick, bins, f_thick_turb, label=L"l_{c} = %$correlation_length", color=:black, alpha=opacity, linestyle=:dash, linewidth=linewidth)
    end

    filename_thick = joinpath(output_dir, "$(turbulence_model)_line_profile_thick_incl$(incl_angle).pdf")
    savefig(plt_thick, filename_thick)
    display(plt_thick)
end