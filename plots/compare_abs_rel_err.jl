using CairoMakie
using CSV, DataFrames
using LaTeXStrings

include("utils.jl")

df_uniform_err = CSV.read(joinpath(@__DIR__, "../pdmk/accuracy/uniform_sys.csv"), DataFrame)
df_nonuniform_err = CSV.read(joinpath(@__DIR__, "../pdmk/accuracy/nonuniform_sys.csv"), DataFrame)

function plot_accuracy!(df, ax, abs_or_rel)
    Ns = unique(df.N)
    err_3_mean, err_6_mean, err_9_mean = zeros(length(Ns)), zeros(length(Ns)), zeros(length(Ns))
    err_3_max, err_6_max, err_9_max = zeros(length(Ns)), zeros(length(Ns)), zeros(length(Ns))
    err_3_min, err_6_min, err_9_min = zeros(length(Ns)), zeros(length(Ns)), zeros(length(Ns))
    err_3_std, err_6_std, err_9_std = zeros(length(Ns)), zeros(length(Ns)), zeros(length(Ns))
    for i in 1:length(Ns)
        N = Ns[i]

        if abs_or_rel == "abs"
            err_3 = df[df.N .== N, :abs_err_3]
            err_6 = df[df.N .== N, :abs_err_6]
            err_9 = df[df.N .== N, :abs_err_9]
        else
            err_3 = df[df.N .== N, :rel_err_3]
            err_6 = df[df.N .== N, :rel_err_6]
            err_9 = df[df.N .== N, :rel_err_9]
        end

        err_3_mean[i] = geometric_mean(err_3)
        err_6_mean[i] = geometric_mean(err_6)
        err_9_mean[i] = geometric_mean(err_9)
        err_3_max[i] = maximum(err_3)
        err_6_max[i] = maximum(err_6)
        err_9_max[i] = maximum(err_9)
        err_3_min[i] = minimum(err_3)
        err_6_min[i] = minimum(err_6)
        err_9_min[i] = minimum(err_9)
        err_3_std[i] = std(err_3)
        err_6_std[i] = std(err_6)
        err_9_std[i] = std(err_9)
    end

    scatter!(ax, Ns, err_3_mean, label="3 digits", color=colors[1], markersize=markersize, marker=markerstyle[1], strokewidth=strokewidth,)
    errorbars!(ax, Ns, err_3_mean, err_3_mean .- err_3_min, err_3_max .- err_3_mean, color = colors[1], whiskerwidth = 10)

    scatter!(ax, Ns, err_6_mean, label="6 digits", color=colors[2], markersize=markersize, marker=markerstyle[2], strokewidth=strokewidth)
    errorbars!(ax, Ns, err_6_mean, err_6_mean .- err_6_min, err_6_max .- err_6_mean, color = colors[2], whiskerwidth = 10)

    scatter!(ax, Ns, err_9_mean, label="9 digits", color=colors[3], markersize=markersize, marker=markerstyle[3], strokewidth=strokewidth)
    errorbars!(ax, Ns, err_9_mean, err_9_mean .- err_9_min, err_9_max .- err_9_mean, color = colors[3], whiskerwidth = 10)

    # axislegend(ax_1, position=:lt)
    ylims!(ax, 10^(-14), 10^(0))
    xlims!(ax, 10^(2.9), 10^(6.1))
end

begin
    fig = Figure(size = (1000, 400), fontsize = 20)
    ax_1 = Axis(fig[1, 1], xlabel = L"N", ylabel = L"\mathcal{E}_r", xscale = log10, yscale = log10)
    ax_2 = Axis(fig[1, 2], xlabel = L"N", ylabel = L"\mathcal{E}", xscale = log10, yscale = log10)

    plot_accuracy!(df_uniform_err, ax_1, "rel")
    plot_accuracy!(df_uniform_err, ax_2, "abs")

    axislegend(ax_1, position = :lb, nbanks = 1, orientation = :horizontal)

    save(joinpath(@__DIR__, "../figs/compare_abs_rel_err_uniform.svg"), fig)
    save(joinpath(@__DIR__, "../figs/compare_abs_rel_err_uniform.png"), fig)

    fig
end