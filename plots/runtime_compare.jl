using CairoMakie
using CSV
using DataFrames
using LaTeXStrings

include("utils.jl")


fmmmc_runtime = CSV.read("../fmmmc/runtime/fmmmc_pbc_local_uniform.csv", DataFrame)
pdmk_runtime  = CSV.read("../pdmk/runtime/pdmk_runtime_uniform_float.csv", DataFrame)

fmmmc_runtime_nu = CSV.read("../fmmmc/runtime/fmmmc_pbc_local_nonuniform.csv", DataFrame)
pdmk_runtime_nu  = CSV.read("../pdmk/runtime/pdmk_runtime_nonuniform_float.csv", DataFrame)

fmmmc_ns = unique(fmmmc_runtime.N)
pdmk_ns  = unique(pdmk_runtime.N)

fmmmc_ns_nu = unique(fmmmc_runtime_nu.N)
pdmk_ns_nu  = unique(pdmk_runtime_nu.N)


begin
    fig = Figure(size=(1000, 450), fontsize=20)
    
  
    ax_args = (
        xlabel = L"N", 
        ylabel = L"$T$ (ms)", 
        xscale = log10, 
        yscale = log10,
        xminorticksvisible = true, xminorgridvisible = true, xminorticks = IntervalsBetween(5),
        yminorticksvisible = true, yminorgridvisible = true, yminorticks = IntervalsBetween(5),
        xticks = (10.0 .^ (4:7), [L"10^4", L"10^5", L"10^6", L"10^7"]),
        yticks = (10.0 .^ (-2.5:0.5:0.5), [L"10^{-2.5}", L"10^{-2.0}", L"10^{-1.5}", L"10^{-1.0}", L"10^{-0.5}", L"10^{0.0}", L"10^{0.5}"])
    )

    ax1 = Axis(fig[1, 1]; ax_args...)
    ax2 = Axis(fig[1, 2]; ax_args...)


    scatter!(ax1, pdmk_ns, pdmk_runtime.time_propose .* 1000, label="DMK-MC propose", color=colors[1], markersize=markersize, marker=markerstyle[1], strokewidth=strokewidth)
    scatter!(ax1, pdmk_ns, pdmk_runtime.time_accept .* 1000, label="DMK-MC accept", color=colors[2], markersize=markersize, marker=markerstyle[2], strokewidth=strokewidth)
    scatter!(ax1, fmmmc_ns, fmmmc_runtime.local_propose .* 1000, label="FMM Local propose", color=colors[3], markersize=markersize, marker=markerstyle[3], strokewidth=strokewidth)
    scatter!(ax1, fmmmc_ns, fmmmc_runtime.local_accept .* 1000, label="FMM Local accept", color=colors[4], markersize=markersize, marker=markerstyle[4], strokewidth=strokewidth)

    xlims!(ax1, 10^(3.9), 10^(7.1))
    ylims!(ax1, 10^(-2.5), 10^(0.5))


    scatter!(ax2, pdmk_ns_nu, pdmk_runtime_nu.time_propose .* 1000, label="DMK-MC propose", color=colors[1], markersize=markersize, marker=markerstyle[1], strokewidth=strokewidth)
    scatter!(ax2, pdmk_ns_nu, pdmk_runtime_nu.time_accept .* 1000, label="DMK-MC accept", color=colors[2], markersize=markersize, marker=markerstyle[2], strokewidth=strokewidth)
    scatter!(ax2, fmmmc_ns_nu, fmmmc_runtime_nu.local_propose .* 1000, label="FMM Local propose", color=colors[3], markersize=markersize, marker=markerstyle[3], strokewidth=strokewidth)
    scatter!(ax2, fmmmc_ns_nu, fmmmc_runtime_nu.local_accept .* 1000, label="FMM Local accept", color=colors[4], markersize=markersize, marker=markerstyle[4], strokewidth=strokewidth)
    
    xlims!(ax2, 10^(3.9), 10^(7.1))
    ylims!(ax2, 10^(-2.5), 10^(0.5))


    Legend(fig[0, :], ax1, orientation=:horizontal, nbanks=1, labelsize = 15)

    text!(ax1, 0, 1, space = :relative, text = "(a)", fontsize = 25, align = (:left, :top), offset = (12, -2), font = :bold)
    text!(ax2, 0, 1, space = :relative, text = "(b)", fontsize = 25, align = (:left, :top), offset = (12, -2), font = :bold)

    save("../figs/runtime_compare_terms_log_log.svg", fig)
    save("../figs/runtime_compare_terms_log_log.pdf", fig)
    fig
end


begin
    fig = Figure(size=(1000, 450), fontsize=20)
    
    ax_args = (
        xlabel = L"N", 
        ylabel = L"$T$ (ms)", 
        xscale = log10, 
        xminorticksvisible = true, xminorgridvisible = true, xminorticks = IntervalsBetween(5),
        yminorticksvisible = true, yminorgridvisible = true, yminorticks = IntervalsBetween(5),
        xticks = (10.0 .^ (4:7), [L"10^4", L"10^5", L"10^6", L"10^7"])
    )

    ax1 = Axis(fig[1, 1]; ax_args...)
    ax2 = Axis(fig[1, 2]; ax_args...)

    scatter!(ax1, pdmk_ns, pdmk_runtime.time_propose .* 1000, label="DMK-MC propose", color=colors[1], markersize=markersize, marker=markerstyle[1], strokewidth=strokewidth)
    scatter!(ax1, pdmk_ns, pdmk_runtime.time_accept .* 1000, label="DMK-MC accept", color=colors[2], markersize=markersize, marker=markerstyle[2], strokewidth=strokewidth)
    scatter!(ax1, fmmmc_ns, fmmmc_runtime.local_propose .* 1000, label="FMM Local propose", color=colors[3], markersize=markersize, marker=markerstyle[3], strokewidth=strokewidth)
    scatter!(ax1, fmmmc_ns, fmmmc_runtime.local_accept .* 1000, label="FMM Local accept", color=colors[4], markersize=markersize, marker=markerstyle[4], strokewidth=strokewidth)

    xlims!(ax1, 10^(3.9), 10^(7.1))
    ylims!(ax1, 0.0, 1.0)

    scatter!(ax2, pdmk_ns_nu, pdmk_runtime_nu.time_propose .* 1000, label="DMK-MC propose", color=colors[1], markersize=markersize, marker=markerstyle[1], strokewidth=strokewidth)
    scatter!(ax2, pdmk_ns_nu, pdmk_runtime_nu.time_accept .* 1000, label="DMK-MC accept", color=colors[2], markersize=markersize, marker=markerstyle[2], strokewidth=strokewidth)
    scatter!(ax2, fmmmc_ns_nu, fmmmc_runtime_nu.local_propose .* 1000, label="FMM Local propose", color=colors[3], markersize=markersize, marker=markerstyle[3], strokewidth=strokewidth)
    scatter!(ax2, fmmmc_ns_nu, fmmmc_runtime_nu.local_accept .* 1000, label="FMM Local accept", color=colors[4], markersize=markersize, marker=markerstyle[4], strokewidth=strokewidth)
    
    xlims!(ax2, 10^(3.9), 10^(7.1))
    ylims!(ax2, 0.0, 2.0)

    Legend(fig[0, :], ax1, orientation=:horizontal, nbanks=1, labelsize = 15)

    text!(ax1, 0, 1, space = :relative, text = "(a)", fontsize = 25, align = (:left, :top), offset = (12, -2), font = :bold)
    text!(ax2, 0, 1, space = :relative, text = "(b)", fontsize = 25, align = (:left, :top), offset = (12, -2), font = :bold)

    save("../figs/runtime_compare_terms_log_y.svg", fig)
    save("../figs/runtime_compare_terms_log_y.pdf", fig)
    fig
end


begin
    fig = Figure(size=(500, 450), fontsize=20)
    ax = Axis(fig[1, 1], xlabel=L"N", ylabel=L"$T$ (ms)", xscale = log10, yticks = 0.0:0.1:1.2) # x-legend 在此

    scatter!(ax, pdmk_ns, pdmk_runtime.time_propose .* 1000, label="DMK-MC propose", color=colors[1], markersize=markersize, marker=markerstyle[1], strokewidth=strokewidth)
    scatter!(ax, pdmk_ns, pdmk_runtime.time_accept .* 1000, label="DMK-MC accept", color=colors[2], markersize=markersize, marker=markerstyle[2], strokewidth=strokewidth)
    scatter!(ax, fmmmc_ns, fmmmc_runtime.local_propose .* 1000, label="FMM Local propose", color=colors[3], markersize=markersize, marker=markerstyle[3], strokewidth=strokewidth)
    scatter!(ax, fmmmc_ns, fmmmc_runtime.local_accept .* 1000, label="FMM Local accept", color=colors[4], markersize=markersize, marker=markerstyle[4], strokewidth=strokewidth)

    xlims!(ax, 10^(3.9), 10^(7.1))
    ylims!(ax, 0.0, 1.0)

    Legend(fig[0, 1], ax, orientation=:horizontal, nbanks=2, labelsize = 15)

    save("../figs/runtime_compare_terms.svg", fig)
    save("../figs/runtime_compare_terms.pdf", fig)
    fig
end


begin
    fig = Figure(size=(500, 450), fontsize=20)
    ax = Axis(fig[1, 1], xlabel=L"N", ylabel=L"$T$ (ms)", xscale = log10, yticks = 0.0:0.1:1.2) # x-legend 在此

    scatter!(ax, pdmk_ns, (pdmk_runtime.time_propose .+ pdmk_runtime.time_accept) .* 1000, label="DMK-MC", color=colors[1], markersize=markersize, marker=markerstyle[1], strokewidth=strokewidth)
    scatter!(ax, fmmmc_ns, (fmmmc_runtime.local_propose .+ fmmmc_runtime.local_accept) .* 1000, label="FMM Local", color=colors[2], markersize=markersize, marker=markerstyle[2], strokewidth=strokewidth)

    xlims!(ax, 10^(3.9), 10^(7.1))
    ylims!(ax, 0.0, 1.1)

    Legend(fig[0, 1], ax, orientation=:horizontal, nbanks=1, labelsize = 15)

    save("../figs/runtime_compare_total.svg", fig)
    save("../figs/runtime_compare_total.pdf", fig)
    fig
end