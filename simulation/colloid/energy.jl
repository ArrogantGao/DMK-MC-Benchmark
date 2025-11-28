using JLD2, Molly
using CSV, DataFrames
using Plots

df_mc = CSV.read(joinpath(@__DIR__, "data/energy_mc_accuracy.csv"), DataFrame)
sys_md = load(joinpath(@__DIR__, "data/md.jld2"))["sys"]

energy_mc = df_mc.E_total
energy_md = [sys_md.loggers.energy.history[i].val for i in 1:length(sys_md.loggers.energy.history)]

fig = plot(energy_mc[2:10:end], label="MC", xlabel="Step", ylabel="Energy", title="Energy Comparison")
plot!(energy_md, label="MD")
ylims!(-1.5 * 1e5, 100.0)