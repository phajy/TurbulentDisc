using Gradus, Plots, Measures, LaTeXStrings

include("updated-velocity-functions.jl") 
include("updated-sound-speed.jl")       

plot_font = "Computer Modern"
default(fontfamily=plot_font,
        linewidth=0.8, label=nothing, 
        grid=false)

M = 1.0                  # Black hole mass
a = 0.998                # Black hole spin parameter
m = KerrMetric(M, a)     # Kerr metric
ratio_thin = 0.3         # Eddington ratio for thin disc
ratio_thick = 0.5        # Eddington ratio for thick disc
h = 10.0                 # Lamppost corona height
correlation_length = 1.0 
mach = 1     
turbulence_model = "fbm" # 'perlin' or 'fbm'        

output_dir = "C:\\Users\\Kate\\project\\TurbulentDisc\\updated_klb_plots\\redshift_maps\\turbulent"
mkpath(output_dir)

function turbulent_redshift(m, x_obs, vel_func, a, M, ratio, alpha, correlation_length, mach)
    g_obs = Gradus.metric(m, x_obs)  
    v_obs = SVector{4, eltype(x_obs)}(1, 0, 0, 0) 

    function _internal_turbulent_redshift(m::AbstractMetric, gp, t)
        v_disc = vel_func(m, gp.x[2], gp.x[4], a, M, ratio, alpha, correlation_length, mach)
        g = Gradus.metric(m, gp.x)
        Gradus.RedshiftFunctions._redshift_dotproduct(g, v_disc, g_obs, v_obs, gp)
    end

    return PointFunction(_internal_turbulent_redshift)
end

function velocity_wrapper(m, r, phis, a, M, ratio, alpha, correlation_length, mach)
    if turbulence_model == "perlin"
        return turbulence_perlin(m, r, phis, a, M, ratio, alpha, correlation_length, mach)
    elseif turbulence_model == "fbm"
        return turbulence_fbm(m, r, phis, a, M, ratio, alpha, correlation_length, mach)
    end
end


function height_profile_exponential(ρ)
    if ρ <= Gradus.isco(m) || ρ > 50.0  # Outer radius = 50M
        return 0  
    end
    return 0.3 * ρ * exp(-ρ / 50.0)  # Scale height proportional to radius
end

for incl_angle in [30, 60, 75]
    println("Generating redshift maps for i=$incl_angle, Mach=$mach...")

    x = SVector(0.0, 1000.0, deg2rad(incl_angle), 0.0)

    # --- THIN DISC ---
    ssd = ShakuraSunyaev(m, eddington_ratio=ratio_thin)  
    pf_thin = turbulent_redshift(m, x, velocity_wrapper, 1.0, 1.0, ratio_thin, 0.1, correlation_length, mach)
    pf_thin = pf_thin ∘ ConstPointFunctions.filter_intersected()

    α, β, img_thin = rendergeodesics(
        m, x, ssd, 2000.0,
        αlims=(-25, 25), βlims=(-20, 20),
        image_width=800, image_height=400,
        verbose=true, pf=pf_thin
    )

    plt_thin = heatmap(α, β, img_thin, aspect_ratio=1,
                       xlabel=L"\alpha", ylabel=L"\beta",
                       colorbar_title="Redshift",
                       colorbar_titleposition=:top,
                       color=:magma)

    filename_thin = joinpath(output_dir, "$(turbulence_model)_turbulent_redshift_thin_incl$(incl_angle)_mach$(mach).pdf")
    savefig(plt_thin, filename_thin)
    println("Saved thin disc redshift map to: $filename_thin")
    display(plt_thin)

    # --- THICK DISC ---
    thick_disc = ThickDisc(height_profile_exponential)
    pf_thick = turbulent_redshift(m, x, velocity_wrapper, 1.0, 1.0, ratio_thick, 0.1, correlation_length, mach)
    pf_thick = pf_thick ∘ ConstPointFunctions.filter_intersected()

    α, β, img_thick = rendergeodesics(
        m, x, thick_disc, 2000.0,
        αlims=(-25, 25), βlims=(-20, 20),
        image_width=800, image_height=400,
        verbose=true, pf=pf_thick
    )

    plt_thick = heatmap(α, β, img_thick, aspect_ratio=1,
                        xlabel=L"\alpha", ylabel=L"\beta",
                        colorbar_title="Redshift",
                        colorbar_titleposition=:top,
                        color=:magma)

    filename_thick = joinpath(output_dir, "$(turbulence_model)_turbulent_redshift_thick_incl$(incl_angle)_mach$(mach).pdf")
    savefig(plt_thick, filename_thick)
    println("Saved thick disc redshift map to: $filename_thick")
    display(plt_thick)
end
