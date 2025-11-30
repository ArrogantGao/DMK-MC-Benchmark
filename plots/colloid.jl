using CSV, DataFrames, CairoMakie, LaTeXStrings

df_mc_na = CSV.read(joinpath(@__DIR__, "../simulation/colloid/data/density_mc_Na.csv"), DataFrame)
df_mc_cl = CSV.read(joinpath(@__DIR__, "../simulation/colloid/data/density_mc_Cl.csv"), DataFrame)
df_md_na = CSV.read(joinpath(@__DIR__, "../simulation/colloid/data/density_md_Na.csv"), DataFrame)
df_md_cl = CSV.read(joinpath(@__DIR__, "../simulation/colloid/data/density_md_Cl.csv"), DataFrame)

df_accuracy = CSV.read(joinpath(@__DIR__, "../simulation/colloid/data/dE_mc_accuracy.csv"), DataFrame)
dE_pme = df_accuracy.dE_pme
dE_pdmk = df_accuracy.dE_pdmk

abs_error = log10.(abs.(dE_pme - dE_pdmk))
rel_error = log10.(abs.(dE_pme - dE_pdmk) ./ abs.(df_accuracy.E_pme / 1000))


begin
    fig = Figure(size=(1000, 400), fontsize=20)
    ax = Axis(fig[1, 1], xlabel=L"r \text{ (nm)}", ylabel=L"\rho \text{ (nm}^{-3}\text{)}", yscale = log10, xminorticksvisible = true, xminorgridvisible = true, xminorticks = IntervalsBetween(5), yminorticksvisible = true, yminorgridvisible = true, yminorticks = IntervalsBetween(5), yticks = (10.0 .^ (-2:1:3), [L"10^{-2}", L"10^{-1}", L"10^{0}", L"10^{1}", L"10^{2}", L"10^{3}"]), xticks = (0.0:0.5:3.0, [L"0.0", L"0.5", L"1.0", L"1.5", L"2.0", L"2.5", L"3.0"]))

    lines!(ax, df_mc_cl.r, df_mc_cl.rho, label="MC Cl", color=:blue, linewidth=4)
    lines!(ax, df_md_cl.r, df_md_cl.rho, label="MD Cl", color=:red, linewidth=4, linestyle = :dash)

    lines!(ax, df_mc_na.r, df_mc_na.rho, label="MC Na", color=:green, linewidth=4)
    lines!(ax, df_md_na.r, df_md_na.rho, label="MD Na", color=:purple, linewidth=4, linestyle = :dash)
    axislegend(ax, position=:rt, nbanks=2)

    xlims!(ax, 0.0, 3.0)
    ylims!(ax, 10^(-2), 10^3)

    ax_2 = Axis(fig[1, 2], xlabel = L"\mathcal{E}_r", ylabel = L"P", yticks = (0:0.2:1, [L"0.0", L"0.2", L"0.4", L"0.6", L"0.8", L"1.0"]), xticks = (-6:1:0, [L"10^{-6}", L"10^{-5}", L"10^{-4}", L"10^{-3}", L"10^{-2}", L"10^{-1}", L"10^0"]), xminorticksvisible = true, xminorgridvisible = true, xminorticks = vcat([i .+ log10.(2:2:9) for i in -6:-1]...), yminorticksvisible = true, yminorgridvisible = true, yminorticks = IntervalsBetween(5))
    
    hist!(ax_2, rel_error, bins=50, normalization = :pdf, color = :blue, strokewidth = 1, strokecolor = :black, label = "3 digits")
    axislegend(ax_2, position=:rt)
    xlims!(ax_2, -6, 0)
    ylims!(ax_2, 0, 1)


    text!(ax, 0, 1, space = :relative, text = "(a)", fontsize = 25, align = (:left, :top), offset = (12, -2), font = :bold)
    text!(ax_2, 0, 1, space = :relative, text = "(b)", fontsize = 25, align = (:left, :top), offset = (12, -2), font = :bold)

    # ylims!(ax, 10^(-3.1), 10^(0.1))
    save(joinpath(@__DIR__, "../figs/colloid_sim.svg"), fig)
    fig
end