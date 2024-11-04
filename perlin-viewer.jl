using CoherentNoise
using Plots

perlin = perlin_2d(;seed=1)

f = inv(0.01)
x = f .* collect(range(0, 50, 100)) ./ 50
y = f .* collect(range(0, 50, 100)) ./ 50

image = [sample(perlin, X, Y) for X in x, Y in y]

heatmap(image)