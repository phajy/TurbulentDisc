using CoherentNoise, Chain, Plots, Colors, Random

function generate_perlin_noise(seed)
    Random.seed!(seed) 

    img = @chain mix(opensimplex2_2d(seed=seed), opensimplex2_2d(seed=seed+1), perlin_2d(seed=seed)) begin
        scale(0.75) 
        gen_image
    end

    # Convert RGB{Float64} matrix to grayscale intensity matrix
    return [Gray(float(c.r)) for c in img] 
end


seed_values = collect(1:16)

# Generate noise maps using different seeds
perlin_noise_maps = [generate_perlin_noise(s) for s in seed_values]

p = plot(
    [heatmap(perlin_noise_maps[i], title="", axis=nothing, frame=:none) for i in 1:16]...,
    layout=(4,4),
    size=(1400, 1400), 
    grid=false
)

output_dir = "C:\\Users\\Kate\\project\\TurbulentDisc\\klb_plots\\noise"
mkpath(output_dir)
output_path = joinpath(output_dir, "perlin.pdf")  
savefig(p, output_path)  

display(p) 
