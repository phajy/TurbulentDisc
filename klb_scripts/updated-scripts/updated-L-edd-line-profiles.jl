using Plots, Gradus, Measures, LaTeXStrings
include("updated-velocity-functions.jl")
include("updated-sound-speed.jl")

plot_font = "Computer Modern"
default(fontfamily=plot_font,
        linewidth=0.8, framestyle=:box, label=nothing, grid=false)
scalefontsizes(1)

a = 0.998
M = 1.0
h = 10.0
outer_radius = 50.0
bins = collect(range(0.1, 2.0, 200))

inclinations = [30, 60, 75]
thin_ratios = [0.1, 0.2, 0.3]
thick_ratios = [0.5, 0.7, 1.0]  
opacity_values = [1.0, 0.7, 0.4]

output_dir = "C:\\Users\\Kate\\project\\TurbulentDisc\\updated_klb_plots\\line_profiles\\eddington"
mkpath(output_dir)

# Thin discs
for incl_angle in inclinations
    x = SVector(0.0, 1_000.0, deg2rad(incl_angle), 0.0)

    plt_thin = plot(
        xlabel = L"\nu/\nu_e",
        ylabel = L"\textrm{Flux \ (Arbitrary \ Units)}",
        legend = :topright,
        left_margin = [5mm 0mm], right_margin = [5mm 0mm],
        top_margin = [5mm 0mm], bottom_margin = [5mm 0mm],
        xlims = (0, 2.0)
    )

    for (idx, ratio) in enumerate(thin_ratios)
        m = KerrMetric(M, a)
        d_thin = ShakuraSunyaev(m, eddington_ratio=ratio)
        ε_thin = emissivity_profile(m, d_thin, LampPostModel(h=h))
        f_thin = lineprofile(m, x, d_thin, ε_thin; bins, method=BinningMethod(), callback=domain_upper_hemisphere(), verbose=true)[2]
        f_thin[end] = 0
        plot!(plt_thin, bins, f_thin, linestyle=:solid, color=:black, alpha=opacity_values[idx], label=L"L/L_{Edd}= %$ratio")
    end

    filename = joinpath(output_dir, "thin_disc_incl$(incl_angle).pdf")
    savefig(plt_thin, filename)
    display(plt_thin)
    println("Saved figure to: $filename")
#end

# Thick discs (using h_scale as proxy for L/L_Edd)
fixed_h_scale = 0.3  # 'base' scale-height factor

for incl_angle in inclinations
    x = SVector(0.0, 1_000.0, deg2rad(incl_angle), 0.0)

    plt_thick = plot(
        xlabel = L"\nu/\nu_e",
        ylabel = L"\textrm{Flux \ (Arbitrary \ Units)}",
        legend = :topright,
        left_margin = [5mm 0mm], right_margin = [5mm 0mm],
        top_margin = [5mm 0mm], bottom_margin = [5mm 0mm],
        xlims = (0, 2.0)
    )

    for (idx, ratio) in enumerate(thick_ratios)
        m = KerrMetric(M, a)

        # h_scale as a proxy for Eddington ratio
        h_scale = fixed_h_scale * ratio
        h_scale = min(h_scale, 1.0)  # physically realistic cap at H/R = 1

        function height_profile_exponential(ρ)
            return h_scale * ρ * exp(-ρ / outer_radius)
        end

        d_thick = ThickDisc(height_profile_exponential)
        ε_thick = emissivity_profile(m, d_thick, LampPostModel(h=h))
        f_thick = lineprofile(m, x, d_thick, ε_thick; bins, method=BinningMethod(), callback=domain_upper_hemisphere(), verbose=true)[2]
        f_thick[end] = 0

        plot!(plt_thick, bins, f_thick, linestyle=:solid, color=:black,
              alpha=opacity_values[idx], label=L"L/L_{Edd}= %$ratio")
    end

    filename = joinpath(output_dir, "thick_disc_incl$(incl_angle).pdf")
    savefig(plt_thick, filename)
    display(plt_thick)
    println("Saved figure to: $filename")
end
end