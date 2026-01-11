include("utils.jl")

df_total_energy_err = CSV.read(joinpath(@__DIR__, "../pdmk/accuracy/total_energy.csv"), DataFrame)

function plot_total_energy_err!(df, ax)
    Ns = unique(df.N)
    err_3_mean, err_6_mean, err_9_mean = zeros(length(Ns)), zeros(length(Ns)), zeros(length(Ns))
    err_3_max, err_6_max, err_9_max = zeros(length(Ns)), zeros(length(Ns)), zeros(length(Ns))
    err_3_min, err_6_min, err_9_min = zeros(length(Ns)), zeros(length(Ns)), zeros(length(Ns))
    err_3_std, err_6_std, err_9_std = zeros(length(Ns)), zeros(length(Ns)), zeros(length(Ns))
    for i in 1:length(Ns)
        N = Ns[i]
        err_3 = df[df.N .== N, :rel_err_3]
        err_6 = df[df.N .== N, :rel_err_6]
        err_9 = df[df.N .== N, :rel_err_9]
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
    ylims!(ax, 10^(-14), 10^(0))
    xlims!(ax, 10^(2.9), 10^(4.1))
end

begin
    fig = Figure(size = (500, 400), fontsize = 20)
    ax = Axis(fig[1, 1], xlabel = L"N", ylabel = L"\mathcal{E}_r", xscale = log10, yscale = log10, title = "total energy relative error")
    plot_total_energy_err!(df_total_energy_err, ax)

    axislegend(ax, position = :lb, nbanks = 1, orientation = :horizontal)

    save(joinpath(@__DIR__, "../figs/total_energy_err.png"), fig)

    fig
end