# This script plots the sound speed function, and radial inflow velocity function defined by Pariev & Bromley 1998.
# The sound speed is plotted as a dashed line, and the radial inflow velocity is plotted as a solid line.
# Combinations of Eddington ratios and alpha values are looped over those in Pariev & Bromley's Figure 1.

include("updated-sound-speed.jl")

using Plots, LaTeXStrings

plot_font = "Computer Modern"
default(fontfamily=plot_font,
        linewidth=0.8, framestyle=:box, label=nothing, 
        grid=false)

scalefontsizes(1)

M = 1.0
rPos = collect(range(1.0, stop=25.0, length=500))

# Define spin values and corresponding transparency levels
spin_values = [0.0, 0.5, 0.9, 0.99, 0.998]
alpha_values = [0.3, 0.5, 0.7, 0.85, 1.0]

plt = plot(layout=(2,2), size=(850, 850))

annot_x = 22 
annot_x2 = 24.5
annot_y1 = 0.25 
annot_y2 = 0.22 

# Loop over Eddington ratios and alpha values
for (i, eddington_ratio) in enumerate([1, 0.5])
    for (j, alpha) in enumerate([0.1, 0.3])
        
        idx = (i-1) * 2 + j  

        plot!(plt[idx], 
            xlabel="r/M", ylabel="v/c",
            legend=:topright,
            xlims=(1, 25), ylims=(0, 0.4),
            grid=false, minorgrid=false
        )

        for (s_idx, a) in enumerate(spin_values)
            m = KerrMetric(M = M, a = a)

            # Plot sound speed (dashed line)
            plot!(plt[idx], rPos, [SoundSpeed(m, r, a, M, eddington_ratio) for r in rPos], 
                linestyle=:dash, color=:black, alpha=alpha_values[s_idx], label="a/M = $a")

            # Plot radial inflow (solid line)
            plot!(plt[idx], rPos, [RadialSpeed(m, r, a, alpha, M, eddington_ratio) for r in rPos], 
                linestyle=:solid, color=:black, alpha=alpha_values[s_idx])
        end

        annotate!(plt[idx], annot_x, annot_y1, text(L"\alpha = %$alpha", 12, :right))
        annotate!(plt[idx], annot_x2, annot_y2, text(L"L = %$(eddington_ratio) L_{Edd}", 12, :right))
    end
end

output_dir = "C:\\Users\\Kate\\project\\TurbulentDisc\\updated_klb_plots\\soundspeed"
mkpath(output_dir)
filename = joinpath(output_dir, "sound_speed_comparison.pdf")
savefig(plt, filename)

display(plt)
