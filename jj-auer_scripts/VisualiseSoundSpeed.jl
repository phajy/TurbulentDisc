using Plots
include("SoundandRadialSpeed.jl")

# Define all parameters
M = 1.0
alpha = 0.1
ratio = 0.1
a = 0.998

# Collect radial locations
rPos = collect(range(1.0, stop=25.0, length=500))

# Create empty plot
plt = plot(
    xlabel = L"r \ / \ M",
    ylabel = L"c_s \ / \ c",
    legend = :topright,
)

# Define the values to iterate through
values = [0.998, 0.99, 0.9, 0.5, 0.1]

# Collect sound speed values
sound_speeds = [[] for _ in values]
for (i, a) in enumerate(values)
    m = KerrMetric(M = M, a = a)
    sound_speeds[i] = [SoundSpeed(m, r, a, M, ratio) for r in rPos]
end

# Collect radial speed values
radial_speeds = [[] for _ in values]
for (i, a) in enumerate(values)
    m = KerrMetric(M = M, a = a)
    radial_speeds[i] = [RadialSpeed(m, r, a, alpha, M, ratio) for r in rPos]
end


colors = palette(:seaborn_colorblind)

#Add values to plot
plot!(rPos, sound_speeds[1], color = colors[1], label = L"a=0.998")
plot!(rPos, sound_speeds[2], color = colors[2], label = L"a=0.99")
plot!(rPos, sound_speeds[3], color = colors[3], label = L"a=0.9")
plot!(rPos, sound_speeds[4], color = colors[4], label = L"a=0.5")
plot!(rPos, sound_speeds[5], color = colors[5], label = L"a=0.1")

plot!(rPos, radial_speeds[1], linestyle=:dash, color = colors[1])
plot!(rPos, radial_speeds[2], linestyle=:dash, color = colors[2])
plot!(rPos, radial_speeds[3], linestyle=:dash, color = colors[3])
plot!(rPos, radial_speeds[4], linestyle=:dash, color = colors[4])
plot!(rPos, radial_speeds[5], linestyle=:dash, color = colors[5])

#savefig(plt, "Other/Figs/Report/Other/sound_radial_speed_spin.pdf")