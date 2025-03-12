# This script plots the sound speed function, and radial inflow velocity function defined by Pariev & Bromley 1998.
# The sound speed is plotted as a dashed line, and the radial inflow velocity is plotted as a solid line.
# Combinations of Eddington ratios and alpha values are looped over those in Pariev & Bromley's Figure 1.

include("updated-sound-speed.jl")

using Plots

plot_font = "Computer Modern"
default(fontfamily=plot_font,
        linewidth=0.8, framestyle=:box, label=nothing, 
        grid=false)

scalefontsizes(1)

M = 1.0
rPos = collect(range(1.0, stop=25.0, length=500))

# Define spin values and corresponding transparency levels (higher spin = more opaque)
spin_values = [0.0, 0.5, 0.9, 0.99, 0.998]
alpha_values = [0.3, 0.5, 0.7, 0.85, 1.0] 

# Loop over Eddington ratios and alpha values (combinations as of Pariev and Bromley 1998)
for eddington_ratio in [1, 0.5]
    for alpha in [0.1, 0.3]
    
        plt = plot(
            xlabel = "Radius (r/M)",
            ylabel = "v/c",
            legend = :topright,
            grid=false,          
            minorgrid=false      
        )

        for (idx, a) in enumerate(spin_values)
            m = KerrMetric(M = M, a = a)
            
            # Plot sound speed (dashed line)
            plot!(rPos, [SoundSpeed(m, r, a, M, eddington_ratio) for r in rPos], 
                linestyle=:dash, 
                color=:navy, 
                alpha=alpha_values[idx], 
                label="")

            # Plot radial inflow (solid line) 
            plot!(rPos, [RadialSpeed(m, r, a, alpha, M, eddington_ratio) for r in rPos], 
                linestyle=:solid, 
                color=:navy,  
                alpha=alpha_values[idx], 
                label="a/M = $a")
        end

        display(plt)

        filename = "alpha_$(alpha)_ratio_$(eddington_ratio).pdf"
        dir = "C:\\Users\\Kate\\project\\TurbulentDisc\\updated_klb_plots\\soundspeed"
        savefig(plt, joinpath(dir, filename))
    end
end
