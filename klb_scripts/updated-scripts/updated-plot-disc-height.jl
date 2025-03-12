using Plots, Gradus, LaTeXStrings

plot_font = "Computer Modern"
default(fontfamily=plot_font,
        linewidth=0.8, framestyle=:box, label=nothing, grid=false)
scalefontsizes(1)

M = 1.0   # Black hole mass
a = 0.998 # Spin parameter
m = KerrMetric(M, a)

inner_radius = Gradus.isco(m) 
outer_radius = 400.0          
r_values = collect(range(inner_radius, stop=outer_radius, length=500))

# --- Define Thick Disc Height Profile (Exponential) ---
function height_profile_exponential(ρ, h_scale, outer_radius)
    return h_scale * ρ * exp(-ρ / outer_radius)
end

h_scales = [0.3, 0.8, 1.5, 2.0]
opacities = [1.0, 0.65, 0.4, 0.2]  

plt = plot(
    xlabel = "r/M",
    ylabel = "h(r)",
    legend = :topright,
    xlims = (inner_radius, outer_radius)
)

for (idx, h_scale) in enumerate(h_scales)
    h_values = [height_profile_exponential(r, h_scale, outer_radius) for r in r_values]
    plot!(plt, r_values, h_values, label=L"h_{scale} = %$h_scale", linestyle=:solid, color=:black, alpha=opacities[idx])
end

display(plt)

output_dir = "C:\\Users\\Kate\\project\\TurbulentDisc\\updated_klb_plots\\disc_height"
mkpath(output_dir)
filename = joinpath(output_dir, "disc_height_vs_radius.pdf")
#savefig(plt, filename)
#println("Saved figure to: $filename")
