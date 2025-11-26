using JLD2
using Molly
using CSV, DataFrames

include(joinpath(@__DIR__, "../MollyExt/rdf.jl"))

sys_mc = load(joinpath(@__DIR__, "NaCl_mc.jld2"))["sys"]
sys_md = load(joinpath(@__DIR__, "NaCl_md.jld2"))["sys"]

N_samples = 400

rdf_mc_y_nacl = Vector{Vector}(undef, N_samples)
rdf_md_y_nacl = Vector{Vector}(undef, N_samples)
rdf_mc_y_nana = Vector{Vector}(undef, N_samples)
rdf_mc_y_clcl = Vector{Vector}(undef, N_samples)
rdf_md_y_nana = Vector{Vector}(undef, N_samples)
rdf_md_y_clcl = Vector{Vector}(undef, N_samples)

for i in 0:N_samples-1
    @show i
    mc_coords = sys_mc.loggers.coords.history[end - i * 10]
    md_coords = sys_md.loggers.coords.history[end - i]

    mc_Na_coords = [[x.val for x in mc_coords[i]] for i in 1:500]
    mc_Cl_coords = [[x.val for x in mc_coords[i]] for i in 501:1000]
    md_Na_coords = [[x.val for x in md_coords[i]] for i in 1:500]
    md_Cl_coords = [[x.val for x in md_coords[i]] for i in 501:1000]

    _, rdf_mc_y_nacl[i + 1] = rdf_type(mc_Na_coords, mc_Cl_coords, 9.4, 200)
    _, rdf_md_y_nacl[i + 1] = rdf_type(md_Na_coords, md_Cl_coords, 9.4, 200)
    _, rdf_mc_y_nana[i + 1] = rdf_type(mc_Na_coords, mc_Na_coords, 9.4, 200)
    _, rdf_mc_y_clcl[i + 1] = rdf_type(mc_Cl_coords, mc_Cl_coords, 9.4, 200)
    _, rdf_md_y_nana[i + 1] = rdf_type(md_Na_coords, md_Na_coords, 9.4, 200)
    _, rdf_md_y_clcl[i + 1] = rdf_type(md_Cl_coords, md_Cl_coords, 9.4, 200)
end

rdf_mc_y_nacl_avg = sum(rdf_mc_y_nacl) / length(rdf_mc_y_nacl)
rdf_md_y_nacl_avg = sum(rdf_md_y_nacl) / length(rdf_md_y_nacl)
rdf_mc_y_nana_avg = sum(rdf_mc_y_nana) / length(rdf_mc_y_nana)
rdf_mc_y_clcl_avg = sum(rdf_mc_y_clcl) / length(rdf_mc_y_clcl)
rdf_md_y_nana_avg = sum(rdf_md_y_nana) / length(rdf_md_y_nana)
rdf_md_y_clcl_avg = sum(rdf_md_y_clcl) / length(rdf_md_y_clcl)

bins = range(0.0, 9.4, length = length(rdf_mc_y_nacl_avg) + 1)
xs = [(bins[i] + bins[i + 1]) / 2 for i in 1:length(rdf_mc_y_nacl_avg)]

df = joinpath(@__DIR__, "data/rdf.csv")
CSV.write(df, DataFrame(r = xs, rdf_mc_y_nacl = rdf_mc_y_nacl_avg, rdf_md_y_nacl = rdf_md_y_nacl_avg, rdf_mc_y_nana = rdf_mc_y_nana_avg, rdf_mc_y_clcl = rdf_mc_y_clcl_avg, rdf_md_y_nana = rdf_md_y_nana_avg, rdf_md_y_clcl = rdf_md_y_clcl_avg))