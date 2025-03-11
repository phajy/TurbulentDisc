using Plots
include("SoundandRadialSpeed.jl")

M = 1.0
alpha = 0.1
ratio = 1.0

rPos = collect(range(1.0, stop=25.0, length=500))

plt_sound = plot(
    xlabel = "Radius (r/M)",
    ylabel = "Sound Speed (v/c)",
    title = "Sound Speed vs. Radius",
    legend = :topright,
)

for a in [0.0, 0.5, 0.9, 0.99, 0.998]
    m = KerrMetric(M = M, a = a)
    plot!(rPos, [SoundSpeed(m, r, a, M, ratio) for r in rPos], label="a/M = $a")
end

display(plt_sound)


plt_radial = plot(
    xlabel = "Radius (r/M)",
    ylabel = "Radial Inflow Speed (v/c)",
    title = "Radial Inflow Speed vs. Radius",
    legend = :topright,
)

for a in [0.0, 0.5, 0.9, 0.99, 0.998]
    m = KerrMetric(M = M, a = a)
    plot!(rPos, [RadialSpeed(m, r, a, alpha, M, ratio) for r in rPos], label="a/M = $a")
end

display(plt_radial)