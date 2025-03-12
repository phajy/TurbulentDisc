# Script to plot non-turbulent line profiles for thin and thick discs with varying spin (and inclination angle)

using Plots, Gradus, Measures, LaTeXStrings
include("updated-velocity-functions.jl")
include("updated-sound-speed.jl")

plot_font = "Computer Modern"
default(fontfamily=plot_font,
        linewidth=0.8, framestyle=:box, label=nothing, grid=false)
scalefontsizes(1)

M = 1.0                  # Black hole mass
ratio_thin = 0.3         # Eddington ratio for thin disc
ratio_thick = 0.5        # Eddington ratio for thick disc
h_scale = 0.3            # Scale height for thick disc
h = 10.0                 # Lamppost corona height
outer_radius = 50.0      # Outer disc radius for thick disc
bins = collect(range(0.1, 2.0, 200))

spin_values = [0.0, 0.5, 0.9, 0.99, 0.998]
opacity_values = [1.0, 0.8, 0.6, 0.4, 0.2]  
inclination_angles = [30, 60, 75]  

output_dir = "C:\\Users\\Kate\\project\\TurbulentDisc\\updated_klb_plots\\line_profiles\\spin"
mkpath(output_dir)

# Loop over inclination angles
for incl_angle in inclination_angles
    x = SVector(0.0, 1_000.0, deg2rad(incl_angle), 0.0)

    plt_thin = plot(
        xlabel = L"\nu/\nu_e",
        ylabel = L"\textrm{Flux \ (Arbitrary \ Units)}",
        legend = :topleft,
        title = L"\textrm{i = %$incl_angle^\circ}",
        titlefont = 11,
        left_margin = [5mm 0mm], right_margin = [5mm 0mm],
        top_margin = [5mm 0mm], bottom_margin = [5mm 0mm],
        xlims = (0, 2.0)
    )

    for (idx, a) in enumerate(spin_values)
        m = KerrMetric(M, a)

        # Define the thin accretion disc (Shakura-Sunyaev)
        d_thin = ShakuraSunyaev(m, eddington_ratio=ratio_thin)
        ε_thin = emissivity_profile(m, d_thin, LampPostModel(h=h))

        f_thin = lineprofile(m, x, d_thin, ε_thin; bins, method = BinningMethod(), callback = domain_upper_hemisphere(), verbose=true)[2]
        f_thin[end] = 0  

        plot!(plt_thin, bins, f_thin, label="a/M = $a", color=:black, alpha=opacity_values[idx], linestyle=:solid)
    end

    filename_thin = joinpath(output_dir, "line_profile_thin_a_variation_incl$(incl_angle).pdf")
    savefig(plt_thin, filename_thin)
    display(plt_thin)

    plt_thick = plot(
        xlabel = L"\nu/\nu_e",
        ylabel = L"\textrm{Flux \ (Arbitrary \ Units)}",
        legend = :topleft,
        left_margin = [5mm 0mm], right_margin = [5mm 0mm],
        top_margin = [5mm 0mm], bottom_margin = [5mm 0mm],
        xlims = (0, 2.0)
    )

    for (idx, a) in enumerate(spin_values)
        m = KerrMetric(M, a)

        # Define the thick accretion disc
        function height_profile_exponential(ρ)
            return 0.3 * ρ * exp(-ρ / outer_radius) # 0.3 here is h_scale
        end

        d_thick = ThickDisc(height_profile_exponential)
        ε_thick = emissivity_profile(m, d_thick, LampPostModel(h=h))

        f_thick = lineprofile(m, x, d_thick, ε_thick; bins, method = BinningMethod(), callback = domain_upper_hemisphere(), verbose=true)[2]
        f_thick[end] = 0  

        plot!(plt_thick, bins, f_thick, label="a/M = $a", color=:black, alpha=opacity_values[idx], linestyle=:solid)
    end

    filename_thick = joinpath(output_dir, "line_profile_thick_a_variation_incl$(incl_angle).pdf")
    savefig(plt_thick, filename_thick)
    display(plt_thick)
end

