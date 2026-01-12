include("utils.jl")

using LsqFit

function plot_accuracy!(df, ax)
    Ns = unique(df.N)
    err_3_mean, err_6_mean, err_9_mean = zeros(length(Ns)), zeros(length(Ns)), zeros(length(Ns))
    err_3_max, err_6_max, err_9_max = zeros(length(Ns)), zeros(length(Ns)), zeros(length(Ns))
    err_3_min, err_6_min, err_9_min = zeros(length(Ns)), zeros(length(Ns)), zeros(length(Ns))
    err_3_std, err_6_std, err_9_std = zeros(length(Ns)), zeros(length(Ns)), zeros(length(Ns))
    for i in 1:length(Ns)
        N = Ns[i]
        # err_3 = df[df.N .== N, :abs_err_3]
        # err_6 = df[df.N .== N, :abs_err_6]
        # err_9 = df[df.N .== N, :abs_err_9]

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
    # errorbars!(ax, Ns, err_3_mean, err_3_mean .- err_3_min, err_3_max .- err_3_mean, color = colors[1], whiskerwidth = 10)

    scatter!(ax, Ns, err_6_mean, label="6 digits", color=colors[2], markersize=markersize, marker=markerstyle[2], strokewidth=strokewidth)
    # errorbars!(ax, Ns, err_6_mean, err_6_mean .- err_6_min, err_6_max .- err_6_mean, color = colors[2], whiskerwidth = 10)

    scatter!(ax, Ns, err_9_mean, label="9 digits", color=colors[3], markersize=markersize, marker=markerstyle[3], strokewidth=strokewidth)
    # errorbars!(ax, Ns, err_9_mean, err_9_mean .- err_9_min, err_9_max .- err_9_mean, color = colors[3], whiskerwidth = 10)

    # axislegend(ax_1, position=:lt)
    ylims!(ax, 10^(-11), 10^(-1))
    xlims!(ax, 10^(2.9), 10^(6.1))
end

function plot_runtime!(df, ax, ymax)
    Ns = unique(df.N)
    time_eval_3 = df.te_3
    time_eval_6 = df.te_6
    time_eval_9 = df.te_9
    time_update_3 = df.tu_3
    time_update_6 = df.tu_6
    time_update_9 = df.tu_9

    # fit the runtime for evaluation, the runtime should be of y = a * log x
    for (i, data) in enumerate([time_eval_3, time_eval_6, time_eval_9]) 
        model(x, p) = p[1] .* log10.(x) .+ p[2]
        p0 = [1.0, 1.0]
        fit = curve_fit(model, Ns, data .* 1000, p0)
        x_fit = logrange(10^(2.9), 10^(6.1), 100)
        y_fit = model(x_fit, fit.param)
        lines!(ax, x_fit, y_fit, color=colors[i], linewidth=2, linestyle = :dot)
    end

    scatter!(ax, Ns, time_eval_3 .* 1000, color=colors[1], markersize=markersize, marker=markerstyle[1], strokewidth=strokewidth)
    scatter!(ax, Ns, time_eval_6 .* 1000, color=colors[2], markersize=markersize, marker=markerstyle[2], strokewidth=strokewidth)
    scatter!(ax, Ns, time_eval_9 .* 1000, color=colors[3], markersize=markersize, marker=markerstyle[3], strokewidth=strokewidth)
    # scatter!(ax, Ns, time_update_3 .* 1000, color=colors[4], markersize=markersize, marker=markerstyle[4], strokewidth=strokewidth)
    # scatter!(ax, Ns, time_update_6 .* 1000, color=colors[5], markersize=markersize, marker=markerstyle[5], strokewidth=strokewidth)
    # scatter!(ax, Ns, time_update_9 .* 1000, color=colors[6], markersize=markersize, marker=markerstyle[6], strokewidth=strokewidth)

    xlims!(ax, 10^(2.9), 10^(6.1))
    ylims!(ax, 0.0, ymax)
end

begin
    df = CSV.read("../pdmk/accuracy/uniform_sys.csv", DataFrame)
    df_nu = CSV.read("../pdmk/accuracy/nonuniform_sys.csv", DataFrame)

    df_runtime = CSV.read("../pdmk/runtime/uniform_sys.csv", DataFrame)
    df_runtime_nu = CSV.read("../pdmk/runtime/nonuniform_sys.csv", DataFrame)

    fig = Figure(size=(900, 800), fontsize=20)
    ax_1 = Axis(fig[1, 1], xlabel=L"N", ylabel=L"\mathcal{E}_r \left(\Delta U \right)", xscale = log10, yscale = log10, xminorticksvisible = true, xminorgridvisible = true, xminorticks = IntervalsBetween(5), yminorticksvisible = true, yminorgridvisible = true, yminorticks = IntervalsBetween(5), yticks = (10.0 .^ (-12:-1), [L"10^{-12}", L"10^{-11}", L"10^{-10}", L"10^{-9}", L"10^{-8}", L"10^{-7}", L"10^{-6}", L"10^{-5}", L"10^{-4}", L"10^{-3}", L"10^{-2}", L"10^{-1}"]), xticks = (10.0 .^ (3:7), [L"10^3", L"10^4", L"10^5", L"10^6", L"10^7"]))
    ax_2 = Axis(fig[1, 2], xlabel=L"N", ylabel=L"\mathcal{E}_r \left(\Delta U \right)", xscale = log10, yscale = log10, xminorticksvisible = true, xminorgridvisible = true, xminorticks = IntervalsBetween(5), yminorticksvisible = true, yminorgridvisible = true, yminorticks = IntervalsBetween(5), yticks = (10.0 .^ (-12:-1), [L"10^{-12}", L"10^{-11}", L"10^{-10}", L"10^{-9}", L"10^{-8}", L"10^{-7}", L"10^{-6}", L"10^{-5}", L"10^{-4}", L"10^{-3}", L"10^{-2}", L"10^{-1}"]), xticks = (10.0 .^ (3:7), [L"10^3", L"10^4", L"10^5", L"10^6", L"10^7"]))

    plot_accuracy!(df, ax_1)
    plot_accuracy!(df_nu, ax_2)
    Legend(fig[0, :], ax_1, orientation=:horizontal, nbanks=1, labelsize = 18)

    ax_3 = Axis(fig[2, 1], xlabel = L"N", ylabel = L"$T$ (ms)", xscale = log10, xminorticksvisible = true, xminorgridvisible = true, xminorticks = IntervalsBetween(5), xticks = (10.0 .^ (3:7), [L"10^3", L"10^4", L"10^5", L"10^6", L"10^7"]), yticks = (0.0:0.1:0.3, [L"0.0", L"0.1", L"0.2", L"0.3"]), yminorticksvisible = true, yminorgridvisible = true, yminorticks = IntervalsBetween(5))
    ax_4 = Axis(fig[2, 2], xlabel = L"N", ylabel = L"$T$ (ms)", xscale = log10, xminorticksvisible = true, xminorgridvisible = true, xminorticks = IntervalsBetween(5), xticks = (10.0 .^ (3:7), [L"10^3", L"10^4", L"10^5", L"10^6", L"10^7"]), yticks = (0.0:0.3:0.9, [L"0.0", L"0.3", L"0.6", L"0.9"]), yminorticksvisible = true, yminorgridvisible = true, yminorticks = IntervalsBetween(5))
    plot_runtime!(df_runtime, ax_3, 0.3)
    plot_runtime!(df_runtime_nu, ax_4, 0.9)

    text!(ax_1, 0, 1, space = :relative, text = "(a)", fontsize = 25, align = (:left, :top), offset = (12, -2), font = :bold)
    text!(ax_2, 0, 1, space = :relative, text = "(b)", fontsize = 25, align = (:left, :top), offset = (12, -2), font = :bold)
    text!(ax_3, 0, 1, space = :relative, text = "(c)", fontsize = 25, align = (:left, :top), offset = (12, -2), font = :bold)
    text!(ax_4, 0, 1, space = :relative, text = "(d)", fontsize = 25, align = (:left, :top), offset = (12, -2), font = :bold)

    save("../figs/accuracy_runtime.svg", fig)

    fig
end