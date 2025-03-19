"""
This script generates a 3D visualisation of geodesic paths in a Kerr spacetime for a lamppost model.
It computes:
- Geodesic trajectories from a point source located at height `h` above the black hole.
- The effect of frame-dragging and light bending around a rapidly spinning Kerr black hole.
- A 3D plot of geodesic paths and the event horizon.
"""
using Plots, Gradus

plot_font = "Computer Modern"
default(fontfamily=plot_font,
        linewidth=0.8, framestyle=:box, label=nothing, grid=true)
scalefontsizes(1)

m = KerrMetric(a = 0.998)
h = 30.0

model = LampPostModel(h=h)
sols = tracegeodesics(
    m,
    model,
    2000.0,
    n_samples = 64
)

plot_paths_3d(sols, legend=false, extent = 50, t_span = 100.0, 
            xlabel="x", ylabel="y", zlabel="z")
plot_horizon_3d!(m)

output_dir = "C:\\Users\\Kate\\project\\TurbulentDisc\\updated_klb_plots\\lamppost_visualisation"
mkpath(output_dir)
filename = joinpath(output_dir, "lamppost_visual_$(h).pdf")
savefig(filename)
