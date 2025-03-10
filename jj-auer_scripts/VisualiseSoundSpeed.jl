using Plots
include("SoundSpeed.jl")

M = 1.0
alpha = 0.1
lum = 1e46

plt = plot(
    xlabel = "Radius (r/M)",
    ylabel = "Sound Speed (v/c)",
    title = "Sound Speed vs. Radius",
    legend = :topright,
)

rPos = collect(range(1.0, stop=25.0, length=500))

for a in [0.7, 0.85, 0.9, 0.99, 0.998]
    m = KerrMetric(M = M, a = a)
    plot!(rPos, [SoundSpeed(m, r, a, M, lum) for r in rPos], label="a/M = $a")
end

display(plt)