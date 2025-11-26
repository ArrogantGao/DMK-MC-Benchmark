using CSV, DataFrames, CairoMakie, LaTeXStrings

df_csv = CSV.read(joinpath(@__DIR__, "../simulation/NaCl/data/rdf.csv"), DataFrame)
df_accuracy = CSV.read(joinpath(@__DIR__, "../simulation/NaCl/data/dE_mc_accuracy.csv"), DataFrame)

dE_pme = df_accuracy.dE_pme
dE_pdmk = df_accuracy.dE_pdmk

abs_error = log10.(abs.(dE_pme - dE_pdmk))
rel_error = log10.(abs.(dE_pme - dE_pdmk) ./ abs.(df_accuracy.E_pme / 1000))

begin
    fig_rdf = Figure(size = (1000, 400), fontsize = 20)
    ax_rdf = Axis(fig_rdf[1, 1], xlabel = L"$r$ (nm)", ylabel = L"\text{RDF}~\left(\text{nm}^{-3}\right)", xticks = (0:0.5:3, [L"0.0", L"0.5", L"1.0", L"1.5", L"2.0", L"2.5", L"3.0"]), yticks = (0:1:5, [L"0.0", L"1.0", L"2.0", L"3.0", L"4.0", L"5.0"]), xminorticksvisible = true, xminorgridvisible = true, xminorticks = IntervalsBetween(5), yminorticksvisible = true, yminorgridvisible = true, yminorticks = IntervalsBetween(5))
    lines!(ax_rdf, df_csv.r, df_csv.rdf_mc_y_nacl, label="MC Na-Cl", color = :blue, linewidth = 4)
    lines!(ax_rdf, df_csv.r, df_csv.rdf_md_y_nacl, label="MD Na-Cl", color = :red, linewidth = 4, linestyle = :dash)
    lines!(ax_rdf, df_csv.r, df_csv.rdf_mc_y_nana, label="MC Na-Na", color = :green, linewidth = 4)
    lines!(ax_rdf, df_csv.r, df_csv.rdf_md_y_nana, label="MD Na-Na", color = :purple, linewidth = 4, linestyle = :dash)
    # lines!(ax_rdf, df_csv.r, df_csv.rdf_mc_y_clcl, label="MC Cl-Cl", color = :orange, linewidth = 4)
    # lines!(ax_rdf, df_csv.r, df_csv.rdf_md_y_clcl, label="MD Cl-Cl", color = :brown, linewidth = 4, linestyle = :dash)
    axislegend(ax_rdf, position=:rt, nbanks = 2)
    xlims!(ax_rdf, 0.0, 3.0)
    ylims!(ax_rdf, 0.0, 5.0)

    ax_2 = Axis(fig_rdf[1, 2], xlabel = L"\mathcal{E}_r", ylabel = L"P", yticks = (0:0.2:1, [L"0.0", L"0.2", L"0.4", L"0.6", L"0.8", L"1.0"]), xticks = (-6:1:0, [L"10^{-6}", L"10^{-5}", L"10^{-4}", L"10^{-3}", L"10^{-2}", L"10^{-1}", L"10^0"]), xminorticksvisible = true, xminorgridvisible = true, xminorticks = vcat([i .+ log10.(2:2:9) for i in -6:-1]...), yminorticksvisible = true, yminorgridvisible = true, yminorticks = IntervalsBetween(5))
    
    hist!(ax_2, rel_error, bins=50, normalization = :pdf, color = :blue, strokewidth = 1, strokecolor = :black, label = "3 digits")
    axislegend(ax_2, position=:rt)
    xlims!(ax_2, -6, -1)
    ylims!(ax_2, 0, 1)

    text!(ax_rdf, 0, 1, space = :relative, text = "(a)", fontsize = 25, align = (:left, :top), offset = (12, -2), font = :bold)
    text!(ax_2, 0, 1, space = :relative, text = "(b)", fontsize = 25, align = (:left, :top), offset = (12, -2), font = :bold)

    save(joinpath(@__DIR__, "../figs/NaCl_sim.svg"), fig_rdf)

    fig_rdf
end