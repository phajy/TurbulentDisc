# Script to plot both the line profiles for accretion discswith varying parameters different disc types 
# (thin, thick, Shakura-Sunyaev) for no-turbulence, or turbulence (Perlin/fBm), spin, Eddington ratio and
# varying turbulence parameters (Mach number, correlation length)


"""
Parameters:

- `a`: Kerr black hole spin parameter
- `M`: Black hole mass
- `ratio`: Eddington ratio
- `alpha`: Alpha viscosity parameter
- `bins`: Number of bins for the line profile
- `ε`: Emissivity profile
- `correlation_length`: Correlation length for the turbulence
- `mach`: Mach number for the turbulence

"""

using Plots, Gradus, Measures, LaTeXStrings
include("updated-velocity-functions.jl")
include("updated-sound-speed.jl")

plot_font = "Computer Modern"
default(fontfamily=plot_font,
        linewidth=0.8, framestyle=:box, label=nothing, grid=false)
scalefontsizes(1)

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

function velocity_wrapper(m, r, phis, a, M, ratio, alpha, correlation_length, mach)
    return turbulence_perlin(m, r, phis, a, M, ratio, alpha, correlation_length, mach)
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

incl_angle = 60
x = SVector(0.0, 1_000.0, deg2rad(incl_angle), 0.0)

a = 0.998                # Black hole spin parameter
M = 1.0                  # Black hole mass
ratio = 0.3              # Eddington ratio
alpha = 0.1              # Alpha viscosity parameter
h = 10.0                 # Lamppost height
mach = 2                 # Mach number 
correlation_length = 1   # Correlation length

m = KerrMetric(M, a)
inner_radius = Gradus.isco(m)
outer_radius = 400.0
d = ShakuraSunyaev(m, eddington_ratio=ratio)
ε = emissivity_profile(m, d, LampPostModel(h=10.0)) 

bins = collect(range(0.1, 2.0, 200))

f_turb = calculate_turbulent_line_profile(m, x, d, bins, correlation_length, a, M, ratio, alpha, mach, ε)
f_laminar = calculate_zero_turbulence_line_profile(m, x, d, bins, ε)

plt = plot(
    xlabel = L"ν/ν_e",
    ylabel = L"\textrm{Flux \ (Arbitrary Units)}",
    legend = :topleft,
    left_margin = [5mm 0mm],
    right_margin = [5mm 0mm],
    top_margin = [5mm 0mm],
    bottom_margin = [5mm 0mm],
    xlims = (0, 2.0)
)

plot!(plt, bins, f_laminar, label=L"\textrm{Zero \ Turbulence}", color=:black, linestyle=:solid, linewidth=0.8)
plot!(plt, bins, f_turb, label=L"\textrm{Turbulence, Mach = %$mach}", color=:black, linestyle=:dash, linewidth=0.8)

display(plt)

output_dir = "C:\\Users\\Kate\\project\\TurbulentDisc\\updated_klb_plots\\line_profiles"
mkpath(output_dir)
filename = joinpath(output_dir, "line_profile_turb_vs_nonturb_a$(a)_mach$(mach).pdf")
#savefig(plt, filename)
#println("Saved figure to: $filename")
