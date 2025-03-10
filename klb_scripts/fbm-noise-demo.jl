using CoherentNoise, Chain, Plots, Colors, Random 

function generate_noise_image(; octaves=4, frequency=1.0, lacunarity=2.0, persistence=0.5)
    params = Dict(:octaves => octaves, :frequency => frequency, :lacunarity => lacunarity, :persistence => persistence)

    # Generate noise image
    img = @chain opensimplex2_3d(seed=1) begin
        fbm_fractal_3d(; source=_, params...)
        gen_image
    end

    # Convert RGB{Float64} matrix to grayscale intensity matrix
    return [Gray(float(c.r)) for c in img]  # Extract red channel as grayscale
end

# Define parameter variations (4 values each)
octaves_values = [2, 4, 6, 10]
frequency_values = [0.5, 1.0, 2.0, 5.0]
lacunarity_values = [1.5, 2.0, 2.5, 5.0]
persistence_values = [0.3, 0.5, 0.7, 2.0]

# Fix middle values for other parameters
default_octaves = 4
default_frequency = 1.0
default_lacunarity = 2.0
default_persistence = 0.5
fixed_seed = 1   

# Generate noise maps for each row
noise_maps = vcat(
    # Row 1: Vary Octaves 
    [begin
        if o == default_octaves Random.seed!(fixed_seed) end  
        generate_noise_image(octaves=o, frequency=default_frequency, lacunarity=default_lacunarity, persistence=default_persistence)
    end for o in octaves_values],

    # Row 2: Vary Frequency
    [begin
        if f == default_frequency Random.seed!(fixed_seed) end
        generate_noise_image(octaves=default_octaves, frequency=f, lacunarity=default_lacunarity, persistence=default_persistence)
    end for f in frequency_values],

    # Row 3: Vary Lacunarity
    [begin
        if l == default_lacunarity Random.seed!(fixed_seed) end
        generate_noise_image(octaves=default_octaves, frequency=default_frequency, lacunarity=l, persistence=default_persistence)
    end for l in lacunarity_values],

    # Row 4: Vary Persistence
    [begin
        if p == default_persistence Random.seed!(fixed_seed) end
        generate_noise_image(octaves=default_octaves, frequency=default_frequency, lacunarity=default_lacunarity, persistence=p)
    end for p in persistence_values]
)


p = plot(
    [heatmap(noise_maps[i], title="", axis=nothing, frame=:none) for i in 1:16]...,
    layout=(4,4), 
    size=(1400, 1400),  
    grid=false 
)


output_dir = "C:\\Users\\Kate\\project\\TurbulentDisc\\klb_plots\\noise"
mkpath(output_dir)
output_path = joinpath(output_dir, "fbm.pdf")  
savefig(p, output_path)

display(p)  
