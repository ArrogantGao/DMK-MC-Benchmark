include("utils.jl")

fmmmc_runtime = CSV.read("../fmmmc/runtime/fmmmc_pbc_local_tset_combined.csv", DataFrame)
pdmk_runtime = CSV.read("../pdmk/runtime/pdmk_runtime_float.csv", DataFrame)

fmmmc_ns = unique(fmmmc_runtime.N)
pdmk_ns = unique(pdmk_runtime.N)

begin
    fig = Figure(size=(500, 450), fontsize=20)
    ax = Axis(fig[1, 1], xlabel=L"N", ylabel=L"$T$ (ms)", xscale = log10, yscale = log10)

    scatter!(ax, pdmk_ns, pdmk_runtime.time_propose .* 1000, label="PDMK propose", color=colors[1], markersize=markersize, marker=markerstyle[1], strokewidth=strokewidth)
    scatter!(ax, pdmk_ns, pdmk_runtime.time_accept .* 1000, label="PDMK accept", color=colors[2], markersize=markersize, marker=markerstyle[2], strokewidth=strokewidth)
    scatter!(ax, fmmmc_ns, fmmmc_runtime.local_propose .* 1000, label="FMM Local propose", color=colors[3], markersize=markersize, marker=markerstyle[3], strokewidth=strokewidth)
    scatter!(ax, fmmmc_ns, fmmmc_runtime.local_accept .* 1000, label="FMM Local accept", color=colors[4], markersize=markersize, marker=markerstyle[4], strokewidth=strokewidth)
    # scatter!(ax, fmmmc_ns, fmmmc_runtime.multipole_propose .* 1000, label="FMM_M propose", color=colors[5], markersize=markersize, marker=markerstyle[5], strokewidth=strokewidth)
    # scatter!(ax, fmmmc_ns, fmmmc_runtime.multipole_accept .* 1000, label="FMM_M accept", color=colors[6], markersize=markersize, marker=markerstyle[6], strokewidth=strokewidth)

    xlims!(ax, 10^(3.9), 10^(7.1))
    ylims!(ax, 10^(-2.5), 10^(0.5))

    Legend(fig[0, 1], ax, orientation=:horizontal, nbanks=2, labelsize = 15)

    save("../figs/runtime_compare_terms_log_log.svg", fig)
    fig
end


begin
    fig = Figure(size=(500, 450), fontsize=20)
    ax = Axis(fig[1, 1], xlabel=L"N", ylabel=L"$T$ (ms)", xscale = log10, yticks = 0.0:0.1:1.2)

    scatter!(ax, pdmk_ns, pdmk_runtime.time_propose .* 1000, label="PDMK propose", color=colors[1], markersize=markersize, marker=markerstyle[1], strokewidth=strokewidth)
    scatter!(ax, pdmk_ns, pdmk_runtime.time_accept .* 1000, label="PDMK accept", color=colors[2], markersize=markersize, marker=markerstyle[2], strokewidth=strokewidth)
    scatter!(ax, fmmmc_ns, fmmmc_runtime.local_propose .* 1000, label="FMM Local propose", color=colors[3], markersize=markersize, marker=markerstyle[3], strokewidth=strokewidth)
    scatter!(ax, fmmmc_ns, fmmmc_runtime.local_accept .* 1000, label="FMM Local accept", color=colors[4], markersize=markersize, marker=markerstyle[4], strokewidth=strokewidth)
    # scatter!(ax, fmmmc_ns, fmmmc_runtime.multipole_propose .* 1000, label="FMM_M propose", color=colors[5], markersize=markersize, marker=markerstyle[5], strokewidth=strokewidth)
    # scatter!(ax, fmmmc_ns, fmmmc_runtime.multipole_accept .* 1000, label="FMM_M accept", color=colors[6], markersize=markersize, marker=markerstyle[6], strokewidth=strokewidth)

    xlims!(ax, 10^(3.9), 10^(7.1))
    ylims!(ax, 0.0, 1.0)

    Legend(fig[0, 1], ax, orientation=:horizontal, nbanks=2, labelsize = 15)

    save("../figs/runtime_compare_terms.svg", fig)
    fig
end

begin
    fig = Figure(size=(500, 450), fontsize=20)
    ax = Axis(fig[1, 1], xlabel=L"N", ylabel=L"$T$ (ms)", xscale = log10, yticks = 0.0:0.1:1.2)

    scatter!(ax, pdmk_ns, (pdmk_runtime.time_propose .+ pdmk_runtime.time_accept) .* 1000, label="PDMK", color=colors[1], markersize=markersize, marker=markerstyle[1], strokewidth=strokewidth)
    scatter!(ax, fmmmc_ns, (fmmmc_runtime.local_propose .+ fmmmc_runtime.local_accept) .* 1000, label="FMM Local", color=colors[2], markersize=markersize, marker=markerstyle[2], strokewidth=strokewidth)
    # scatter!(ax, fmmmc_ns, (fmmmc_runtime.multipole_propose .+ fmmmc_runtime.multipole_accept) .* 1000, label="FMM Multipole", color=colors[3], markersize=markersize, marker=markerstyle[3], strokewidth=strokewidth)

    xlims!(ax, 10^(3.9), 10^(7.1))
    ylims!(ax, 0.0, 1.1)

    Legend(fig[0, 1], ax, orientation=:horizontal, nbanks=1, labelsize = 15)

    save("../figs/runtime_compare_total.svg", fig)
    fig
end