# This script plots polar heatmap of specific turbulence model, correlation length and Mach number

using Plots, Gradus, Measures, LaTeXStrings

include("updated-velocity-functions.jl") 
include("updated-sound-speed.jl")       

plot_font = "Computer Modern"
default(fontfamily=plot_font,
        linewidth=0.8, framestyle=:box, label=nothing, grid=false)
scalefontsizes(1)

function logrange(first, last, num)
    return 10 .^ collect(range(log10(first), log10(last), length=num))
end

M = 1.0                  # Black hole mass
a = 0.998                # Black hole spin parameter
outer_radius = 50.0      # Outer disc outer radius
h_scale = 0.3            # Scale height for thick disc
bins = 2000              
turbulence_model = "perlin" # 'perlin' or 'fbm'
alpha = 0.1
ratio = 0.3
correlation_length = 1
mach = 5 
m = KerrMetric(M, a)      
inner_radius = Gradus.isco(m)  

radii = logrange(inner_radius, outer_radius, bins)
θ = collect(range(0, 2π, bins))  # Angular grid

function velocity_wrapper(m, r, θ, a, M, correlation_length, mach)
    if turbulence_model == "perlin"
        return turbulence_perlin(m, r, θ, a, M, ratio, alpha, correlation_length, mach)
    elseif turbulence_model == "fbm"
        return turbulence_fbm(m, r, θ, a, M, ratio, alpha, correlation_length, mach)
    else
        error("Invalid turbulence model specified: $turbulence_model")
    end
end

keplerian_velocities = Gradus.CircularOrbits.fourvelocity.(m, radii)

turbulent_velocities = [
    [velocity_wrapper(m, radii[i], θ[j], a, M, correlation_length, mach) for j in 1:length(θ)]
    for i in 1:length(radii)
]

velocity_differences = [
    [turbulent_velocities[i][j] - keplerian_velocities[i] for j in 1:length(θ)]
    for i in 1:length(radii)
]

difference_mags = [
    [sqrt(sum(x -> x^2, velocity_differences[i][j][2:4])) for j in 1:length(θ)]
    for i in 1:length(radii)
]

keplerian_mags = [
    sqrt(sum(x -> x^2, keplerian_velocities[i][2:4])) for i in 1:length(radii)
]

difference_fracs = [
    [difference_mags[i][j] / keplerian_mags[i] for j in 1:length(θ)]
    for i in 1:length(radii)
]

difference_field = [100 * difference_fracs[i][j] for i in 1:length(radii), j in 1:length(θ)]

plt = heatmap(
    θ,
    log10.(radii), 
    difference_field,
    projection=:polar,
    colorbar_title=L"v_{turb} - v_{Keplerian} \, (\%)",
    colorbar_titlefontsize=14,
    labelpad=5,
    color=:magma,
    clims=(0, 50),  
    grid=true,
    yticks=false,
)


output_dir = "C:\\Users\\Kate\\project\\TurbulentDisc\\updated_klb_plots\\heatmaps"
mkpath(output_dir)
filename = "$(turbulence_model)_heatmap_corr$(correlation_length)_mach$(mach).pdf"
save_path = joinpath(output_dir, filename)
savefig(plt, save_path)


