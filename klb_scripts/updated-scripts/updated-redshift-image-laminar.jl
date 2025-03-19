"""
This script generates redshift images for different (laminar) accretion disc geometries (thin and thick discs).
It computes:
- The redshift map of an accretion disc around a Kerr black hole.
- A Shakura-Sunyaev thin disc model with an Eddington ratio of 0.3.
- The observer's view of the disc at a specified inclination angle.
- Geodesic rendering of the redshift distribution across the observed image plane.
"""
using Gradus, Plots, Measures, LaTeXStrings

plot_font = "Computer Modern"
default(fontfamily=plot_font,
        linewidth=0.8, label=nothing, 
        grid=false)

# Define metric and observer
m = KerrMetric(M=1.0, a=0.998)
incl_angle = 75
x = SVector(0.0, 1000.0, deg2rad(incl_angle), 0.0)

# Define Shakura-Sunyaev disc (thin disc) with L_edd = 0.3
ssd = ShakuraSunyaev(m, eddington_ratio=0.3)
pf = ConstPointFunctions.redshift(m, x) ∘ ConstPointFunctions.filter_intersected()

# Render geodesic redshift map
α, β, img = rendergeodesics(
    m, x, ssd, 2000.0,
    αlims=(-25, 25), βlims=(-20, 20),
    image_width=800, image_height=400,
    verbose=true, pf=pf
)

# Plot redshift image
plt = heatmap(α, β, img, aspect_ratio=1,
              xlabel=L"\alpha",
              ylabel=L"\beta",
              colorbar_title="Redshift",
              colorbar_titleposition=:top,
)

output_dir = "C:\\Users\\Kate\\project\\TurbulentDisc\\updated_klb_plots\\redshift_maps"
mkpath(output_dir)
filename = joinpath(output_dir, "ssd_redshift_$(incl_angle).pdf")
savefig(plt, filename)

println("Saved figure to: $filename")

display(plt)