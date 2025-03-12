# Script to compare thin vs thick disc line profiles as a baseline

using Plots, Gradus, Measures, LaTeXStrings
include("updated-velocity-functions.jl")
include("updated-sound-speed.jl")

plot_font = "Computer Modern"
default(fontfamily=plot_font,
        linewidth=0.8, framestyle=:box, label=nothing, grid=false)
scalefontsizes(1)

# Function for Zero-Turbulence (Laminar) Line Profile
function calculate_zero_turbulence_line_profile(m, x, d, bins, ε)
    _, f = lineprofile(m, x, d, ε; bins, method = BinningMethod(), callback = domain_upper_hemisphere(), verbose = true)
    f[end] = 0
    return f
end

incl_angle = 30
x = SVector(0.0, 1_000.0, deg2rad(incl_angle), 0.0)

a = 0.998                # Black hole spin parameter
M = 1.0                  # Black hole mass
ratio = 0.3              # Eddington ratio
alpha = 0.1              # Alpha viscosity parameter
h = 10.0                 # Lamppost height

m = KerrMetric(M, a)
inner_radius = Gradus.isco(m)
outer_radius = 400.0
bins = collect(range(0.1, 2.0, 200))

# ---- Thin Disc (Shakura-Sunyaev) ----
thin_disc = ShakuraSunyaev(m, eddington_ratio=ratio)
ε_thin = emissivity_profile(m, thin_disc, LampPostModel(h=h))
f_thin = calculate_zero_turbulence_line_profile(m, x, thin_disc, bins, ε_thin)

# ---- Thick Disc (Exponential Height Profile) ----
h_scales = [0.3, 0.8, 1.5, 2.0]  # Different thick disc heights
opacities = [1.0, 0.7, 0.5, 0.3] 

function height_profile_exponential(ρ, inner_radius, outer_radius, h_scale)
    if ρ < inner_radius || ρ > outer_radius
        return -1.0  # Outside the disc
    else
        return h_scale * ρ * exp(-ρ / outer_radius)
    end
end

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

plot!(plt, bins, f_thin, label=L"Thin \ Disc \ (Shakura-Sunyaev)", color=:black, linestyle=:solid, linewidth=0.8)

# Loop over thick disc heights
for (idx, h_scale) in enumerate(h_scales)
    thick_disc = ThickDisc() do ρ
        height_profile_exponential(ρ, inner_radius, outer_radius, h_scale)
    end
    ε_thick = emissivity_profile(m, thick_disc, LampPostModel(h=h))
    f_thick = calculate_zero_turbulence_line_profile(m, x, thick_disc, bins, ε_thick)
    
    plot!(plt, bins, f_thick, label=L"Thick \ Disc , \ h_{scale} = %$h_scale", 
      color=:black, linestyle=:dash, linewidth=0.8, alpha=opacities[idx])
end

display(plt)

output_dir = "C:\\Users\\Kate\\project\\TurbulentDisc\\updated_klb_plots\\line_profiles"
mkpath(output_dir)
filename = joinpath(output_dir, "baseline_thin_vs_thick_hscale_variation.pdf")
savefig(plt, filename)
println("Saved figure to: $filename")
