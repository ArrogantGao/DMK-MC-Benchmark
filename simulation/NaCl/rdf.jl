using JLD2
using CairoMakie
using Molly
using CSV, DataFrames

include(joinpath(@__DIR__, "../MollyExt/rdf.jl"))

sys_mc = load(joinpath(@__DIR__, "NaCl_mc.jld2"))["sys"]
sys_md = load(joinpath(@__DIR__, "NaCl_md.jld2"))["sys"]

energy_mc = CSV.read(joinpath(@__DIR__, "energy_mc.csv"), DataFrame)

md_energy = [a.val for a in sys_md.loggers.energy.history]


begin
    fig = Figure()
    ax = Axis(fig[1, 1])
    lines!(ax, 0:10:2000, md_energy, label="MD")
    lines!(ax, energy_mc.E_total, label="MC")
    axislegend(ax, position=:lt)
    fig
end

rdf_mc_y_nacl = Vector{Vector}(undef, 10)
rdf_md_y_nacl = Vector{Vector}(undef, 10)

for i in 1:10
    _, rdf_mc_y_nacl[i] = rdf_type(sys_mc.loggers.coords.history[end - i * 2][1:500], sys_mc.loggers.coords.history[end - i][501:1000], 9.4u"nm", 200)
    _, rdf_md_y_nacl[i] = rdf_type(sys_md.loggers.coords.history[end - i][1:500], sys_md.loggers.coords.history[end - i][501:1000], 9.4u"nm", 200)
end

rdf_mc_y_nacl_avg = sum(rdf_mc_y_nacl) / length(rdf_mc_y_nacl)
rdf_md_y_nacl_avg = sum(rdf_md_y_nacl) / length(rdf_md_y_nacl)


begin
    bins = range(0.0u"nm", 9.4u"nm", length = length(rdf_mc_y_nacl_avg) + 1)
    xs = [(bins[i] + bins[i + 1]) / 2 for i in 1:length(rdf_mc_y_nacl_avg)]
    fig_rdf = Figure()
    ax_rdf = Axis(fig_rdf[1, 1])
    # lines!(ax_rdf, mc_x_na, avg_rdf_mc_y_na, label="MC Na")
    # lines!(ax_rdf, md_x_na, avg_rdf_md_y_na, label="MD Na")
    lines!(ax_rdf, xs, [i.val for i in rdf_mc_y_nacl_avg], label="MC Na-Cl")
    lines!(ax_rdf, xs, [i.val for i in rdf_md_y_nacl_avg], label="MD Na-Cl")
    axislegend(ax_rdf, position=:rt)
    # xlims!(ax_rdf, 0u"nm", 1.0u"nm")
    # ylims!(ax_rdf, 10^(-0.5), 10^(1))
    fig_rdf
end
