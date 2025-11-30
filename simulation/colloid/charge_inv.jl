using Molly, JLD2
using CSV, DataFrames
using Plots

function coords_to_val(coords)
    return [[[x.val for x in c] for c in coord] for coord in coords]
end

include(joinpath(@__DIR__, "../MollyExt/charge_int.jl"))
include(joinpath(@__DIR__, "../MollyExt/rdf.jl"))
sys_mc = load(joinpath(@__DIR__, "data/mc.jld2"))["sys"]
sys_md = load(joinpath(@__DIR__, "data/md.jld2"))["sys"]

coords_mc = sys_mc.loggers.coords.history
coords_md = sys_md.loggers.coords.history

charges_mc = [atom.charge for atom in sys_mc.atoms]
charges_md = [atom.charge for atom in sys_md.atoms]
coords_mc_val = coords_to_val(coords_mc[end - 8000:1:end])
coords_md_val = coords_to_val(coords_md[end - 8000:1:end])


L = 9.4
npts_cl = 500
xbins_cl = range(0.0, L / 2, length=npts_cl + 1)
xs_cl = [(xbins_cl[i] + xbins_cl[i + 1]) / 2 for i in 1:length(xbins_cl) - 1]
_, density_mc_Cl = distance_to_center([c[402:end] for c in coords_mc_val], L, L / 2, npts_cl)
_, density_md_Cl = distance_to_center([c[402:end] for c in coords_md_val], L, L / 2, npts_cl)

fig = plot(xs_cl, log10.(density_mc_Cl), label="MC Cl")
plot!(xs_cl, log10.(density_md_Cl), label="MD Cl")
xlims!(0.0, 3.0)
fig

CSV.write(joinpath(@__DIR__, "data/density_mc_Cl.csv"), DataFrame(r = xs_cl, rho = density_mc_Cl))
CSV.write(joinpath(@__DIR__, "data/density_md_Cl.csv"), DataFrame(r = xs_cl, rho = density_md_Cl))

npts_na = 200
xbins_na = range(0.0, L, length=npts_na + 1)
xs_na = [(xbins_na[i] + xbins_na[i + 1]) / 2 for i in 1:length(xbins_na) - 1]
_, density_mc_Na = distance_to_center([c[2:401] for c in coords_mc_val], L, L, npts_na)
_, density_md_Na = distance_to_center([c[2:401] for c in coords_md_val], L, L, npts_na)

fig = plot(xs_na, (density_mc_Na), label="MC Na")
plot!(xs_na, (density_md_Na), label="MD Na")
xlims!(0.0, 3.0)
fig

CSV.write(joinpath(@__DIR__, "data/density_mc_Na.csv"), DataFrame(r = xs_na, rho = density_mc_Na))
CSV.write(joinpath(@__DIR__, "data/density_md_Na.csv"), DataFrame(r = xs_na, rho = density_md_Na))