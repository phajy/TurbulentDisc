using Gradus, Plots

# metric and metric parameters
m = KerrMetric(M=1.0, a=1.0)
# observer position
x = SVector(0.0, 1000.0, deg2rad(80), 0.0)
# accretion disc
d = ThinDisc(1.0, 50.0)

# define point function which filters geodesics that intersected the accretion disc
# and use those to calculate redshift
pf = ConstPointFunctions.redshift(m, x) ∘ ConstPointFunctions.filter_intersected()

α, β, img = rendergeodesics(
    m,
    x,
    d,
    # maximum integration time
    2000.0,
    αlims = (-60, 60), 
    βlims = (-30, 35),
    image_width = 800,
    image_height = 400,
    verbose = true,
    pf = pf,
)

heatmap(α, β, img, aspect_ratio = 1, color=:bluesreds)