# Script to generate redshift images of accretion disc geometries (thin and thick discs)

using Gradus, Plots, Measures, LaTeXStrings

plot_font = "Computer Modern"
default(fontfamily=plot_font,
        linewidth=0.8, label=nothing, 
        grid=false)

m = KerrMetric(M=1.0, a=0.998)
x = SVector(0.0, 1000.0, deg2rad(30), 0.0)

h_scale = 0.3
outer_radius = 50.0 

# Define Thick Disc Height Profile
function height_profile_exponential(u)
    r = u isa Number ? u : u[2]  
    if r > outer_radius
        return -1.0
    else
        return h_scale * r * exp(-r / outer_radius)  
    end
end

thick_disc = ThickDisc(height_profile_exponential)
pf = ConstPointFunctions.redshift(m, x) ∘ ConstPointFunctions.filter_intersected()

α, β, img = rendergeodesics(
    m, x, thick_disc, 2000.0,
    αlims=(-25, 25), βlims=(-20, 20),
    image_width=800, image_height=400,
    verbose=true, pf=pf
)

plt = heatmap(α, β, img, aspect_ratio=1,
              xlabel=L"\alpha",
              ylabel=L"\beta",
              colorbar_title="Redshift",
              colorbar_titleposition=:top,
)

output_dir = "C:\\Users\\Kate\\project\\TurbulentDisc\\updated_klb_plots\\redshift_maps"
mkpath(output_dir)
filename = joinpath(output_dir, "thick_disc_redshift.pdf")
savefig(plt, filename)

println("Saved figure to: $filename")

display(plt)