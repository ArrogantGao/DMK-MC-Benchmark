include("utils.jl")

df_nacl = CSV.read(joinpath(@__DIR__, "../simulation/NaCl/data/dE_mc_accuracy.csv"), DataFrame)
df_colloid = CSV.read(joinpath(@__DIR__, "../simulation/colloid/data/dE_mc_accuracy.csv"), DataFrame)

# scatter the accept rate vs the relative error of the accept rate

dE_PME_nacl = df_nacl.dE_pme
dE_PDMK_nacl = df_nacl.dE_pdmk
dE_PME_colloid = df_colloid.dE_pme
dE_PDMK_colloid = df_colloid.dE_pdmk

KbT = 2.494338785445972
P_PME_nacl = exp.( - dE_PME_nacl / KbT)
P_PDMK_nacl = exp.( - dE_PDMK_nacl / KbT)
P_PME_colloid = exp.( - dE_PME_colloid / KbT)
P_PDMK_colloid = exp.( - dE_PDMK_colloid / KbT)

rel_error_dE_nacl = (abs.(dE_PDMK_nacl - dE_PME_nacl) ./ abs.(dE_PME_nacl))
rel_error_dE_colloid = (abs.(dE_PDMK_colloid - dE_PME_colloid) ./ abs.(dE_PME_colloid))

abs_error_dE_nacl = abs.(dE_PDMK_nacl - dE_PME_nacl)
abs_error_dE_colloid = abs.(dE_PDMK_colloid - dE_PME_colloid)

rel_error_P_nacl = (abs.(P_PDMK_nacl - P_PME_nacl) ./ abs.(P_PME_nacl))
rel_error_P_colloid = (abs.(P_PDMK_colloid - P_PME_colloid) ./ abs.(P_PME_colloid))

avg_P_nacl = mean(rel_error_P_nacl)
avg_P_colloid = mean(rel_error_P_colloid)

@show avg_P_nacl, avg_P_colloid

begin
    fig = Figure(size=(1000, 400), fontsize=20)
    ax_1 = Axis(fig[1, 1], xlabel = "relative error of P", ylabel = "accept rate", xscale = log10, title = "NaCl")
    ax_2 = Axis(fig[1, 2], xlabel = "relative error of P", ylabel = "accept rate", xscale = log10, title = "colloid")

    scatter!(ax_1, rel_error_P_nacl[dE_PME_nacl .> 0], P_PME_nacl[dE_PME_nacl .> 0], color = :blue)
    scatter!(ax_2, rel_error_P_colloid[dE_PME_colloid .> 0], P_PME_colloid[dE_PME_colloid .> 0], color = :red)

    ylims!(ax_1, -0.1, 1.1)
    ylims!(ax_2, -0.1, 1.1)

    save(joinpath(@__DIR__, "../figs/accept_energy_nacl.png"), fig)

    fig
end

begin
    fig = Figure(size=(1000, 400), fontsize=20)
    ax_1 = Axis(fig[1, 1], xlabel = "absolute error of dE", ylabel = "dE", xscale = log10, title = "NaCl")
    ax_2 = Axis(fig[1, 2], xlabel = "absolute error of dE", ylabel = "dE", xscale = log10, title = "colloid")

    scatter!(ax_1, abs_error_dE_nacl[dE_PME_nacl .> 0], dE_PME_nacl[dE_PME_nacl .> 0], color = :blue)
    scatter!(ax_2, abs_error_dE_colloid[dE_PME_colloid .> 0], dE_PME_colloid[dE_PME_colloid .> 0], color = :red)

    # ylims!(ax_1, -0.1, 1.1)
    # ylims!(ax_2, -0.1, 1.1)

    save(joinpath(@__DIR__, "../figs/accept_energy_nacl.png"), fig)

    fig
end