using Plots

x = range(0, 10, 100)

# y = x.^2
y = sin.(x * 0.5)

plot(x, y)
