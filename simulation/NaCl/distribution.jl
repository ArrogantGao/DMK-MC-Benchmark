using CSV, DataFrames
using CairoMakie
using LaTeXStrings

df = CSV.read(joinpath(@__DIR__, "dE_mc_accuracy.csv"), DataFrame)

dE_pme = df.dE_pme
dE_pdmk = df.dE_pdmk

abs_error = log10.(abs.(dE_pme - dE_pdmk))
rel_error = log10.(abs.(dE_pme - dE_pdmk) ./ abs.(df.E_pme / 1000))

begin
    fig = Figure(size = (1000, 400), fontsize = 20)
    ax_1 = Axis(fig[1, 1], xlabel = L"\mathcal{E}", ylabel = L"P", yticks = (0:0.2:1, [L"0.0", L"0.2", L"0.4", L"0.6", L"0.8", L"1.0"]), xticks = (-6:1:0, [L"10^{-6}", L"10^{-5}", L"10^{-4}", L"10^{-3}", L"10^{-2}", L"10^{-1}", L"10^0"]))
    ax_2 = Axis(fig[1, 2], xlabel = L"\mathcal{E}_r", ylabel = L"P", yticks = (0:0.2:1, [L"0.0", L"0.2", L"0.4", L"0.6", L"0.8", L"1.0"]), xticks = (-6:1:0, [L"10^{-6}", L"10^{-5}", L"10^{-4}", L"10^{-3}", L"10^{-2}", L"10^{-1}", L"10^0"]))
    hist!(ax_1, abs_error, bins=50, label="Absolute Error", normalization = :pdf, color = :blue, strokewidth = 1, strokecolor = :black)
    hist!(ax_2, rel_error, bins=50, label="Relative Error", normalization = :pdf, color = :blue, strokewidth = 1, strokecolor = :black)
    axislegend(ax_1, position=:lt)
    axislegend(ax_2, position=:lt)
    xlims!(ax_1, -6, 0)
    xlims!(ax_2, -6, -1)
    ylims!(ax_1, 0, 1)
    ylims!(ax_2, 0, 1)
    fig
end