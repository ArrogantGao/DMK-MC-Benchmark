using Molly, JLD2
using CSV, DataFrames
using Plots

function coords_to_val(coords)
    return [[[x.val for x in c] for c in coord] for coord in coords]
end

include(joinpath(@__DIR__, "../MollyExt/charge_int.jl"))
sys_mc = load(joinpath(@__DIR__, "data/mc.jld2"))["sys"]
sys_md = load(joinpath(@__DIR__, "data/md.jld2"))["sys"]

coords_mc = sys_mc.loggers.coords.history
coords_md = sys_md.loggers.coords.history

charges_mc = [atom.charge for atom in sys_mc.atoms]
charges_md = [atom.charge for atom in sys_md.atoms]

L = 9.4
bins = 100

coords_mc_val = coords_to_val(coords_mc[end - 5000:1:end])
coords_md_val = coords_to_val(coords_md[end - 500:1:end])

xbins = range(0.0, L, length=bins + 1)
xs = [(xbins[i] + xbins[i + 1]) / 2 for i in 1:length(xbins) - 1]

_, density_mc_Na = distance_to_center([c[2:201] for c in coords_mc_val], L, bins)
_, density_mc_Cl = distance_to_center([c[202:701] for c in coords_mc_val], L, bins)
_, density_md_Na = distance_to_center([c[2:201] for c in coords_md_val], L, bins)
_, density_md_Cl = distance_to_center([c[202:701] for c in coords_md_val], L, bins)

fig = plot(xs, log10.(density_mc_Na), label="MC Na")
plot!(xs, log10.(density_mc_Cl), label="MC Cl")
plot!(xs, log10.(density_md_Na), label="MD Na")
plot!(xs, log10.(density_md_Cl), label="MD Cl")
fig

# savefig(fig, joinpath(@__DIR__, "data/density.png"))