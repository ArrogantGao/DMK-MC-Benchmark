using PDMK4MC
using Random
using CSV, DataFrames

Random.seed!(2345)

function accuracy_total_energy(N, ns, L, df, n_trials)

    for i in 1:n_trials

        coords = rand(Float64, 3, N) .* L
        charges = [(-one(Float64))^i for i in 1:N]

        x = Float64.(coords[1, :])
        y = Float64.(coords[2, :])
        z = Float64.(coords[3, :])
        q = ComplexF64.(charges)

        @info "Creating trees"
        tree_3 = PDMK4MC.create_tree(coords, charges; params=PDMK4MC.HPDMKParams(L = L, n_per_leaf=ns, init = DIRECT, digits = 3))
        @info "Tree 3 created"
        tree_6 = PDMK4MC.create_tree(coords, charges; params=PDMK4MC.HPDMKParams(L = L, n_per_leaf=ns, init = DIRECT, digits = 6))
        @info "Tree 6 created"
        tree_9 = PDMK4MC.create_tree(coords, charges; params=PDMK4MC.HPDMKParams(L = L, n_per_leaf=ns, init = DIRECT, digits = 9))
        @info "Tree 9 created"
        tree_12 = PDMK4MC.create_tree(coords, charges; params=PDMK4MC.HPDMKParams(L = L, n_per_leaf=ns, init = DIRECT, digits = 12))
        @info "Tree 12 created"

        n_levels = tree_depth(tree_3)

        for (i, tree) in enumerate([tree_3, tree_6, tree_9, tree_12])
            @info "Forming outgoing PW $i"
            PDMK4MC.form_outgoing_pw!(tree)
            @info "Outgoing PW $i formed"
            PDMK4MC.form_incoming_pw!(tree)
            @info "Incoming PW $i formed"
        end

        E_dmk_3 = PDMK4MC.eval_energy(tree_3)
        E_dmk_6 = PDMK4MC.eval_energy(tree_6)
        E_dmk_9 = PDMK4MC.eval_energy(tree_9)
        E_dmk_12 = PDMK4MC.eval_energy(tree_12)

        abs_err_3 = abs(E_dmk_3 - E_dmk_12)
        abs_err_6 = abs(E_dmk_6 - E_dmk_12)
        abs_err_9 = abs(E_dmk_9 - E_dmk_12)
        rel_err_3 = abs(E_dmk_3 - E_dmk_12) / abs(E_dmk_12)
        rel_err_6 = abs(E_dmk_6 - E_dmk_12) / abs(E_dmk_12)
        rel_err_9 = abs(E_dmk_9 - E_dmk_12) / abs(E_dmk_12)

        CSV.write(df, DataFrame(N = N, ns = ns, L = L, i = i, abs_err_3 = abs_err_3, rel_err_3 = rel_err_3, abs_err_6 = abs_err_6, rel_err_6 = rel_err_6, abs_err_9 = abs_err_9, rel_err_9 = rel_err_9), append=true)

        @info "Total energy error: $(abs_err_3), $(abs_err_6), $(abs_err_9)"
        @info "Total energy relative error: $(rel_err_3), $(rel_err_6), $(rel_err_9)"

        PDMK4MC.destroy_tree!(tree_3)
        PDMK4MC.destroy_tree!(tree_6)
        PDMK4MC.destroy_tree!(tree_9)
        PDMK4MC.destroy_tree!(tree_12)
        GC.gc()
    end

    return nothing
end

function main()
    df = joinpath(@__DIR__, "total_energy.csv")
    CSV.write(df, DataFrame(N = Int[], ns = Int[], L = Float64[], i = Int[], abs_err_3 = Float64[], rel_err_3 = Float64[], abs_err_6 = Float64[], rel_err_6 = Float64[], abs_err_9 = Float64[], rel_err_9 = Float64[]))
    
    rho = 1e4 / 64.63304^3

    accuracy_total_energy(1000, 50, (1000/rho)^(1/3), df, 20)

    N0 = 2000
    for i in 0:3
        N = N0 * 2^i
        L = (N / rho)^(1/3)
        accuracy_total_energy(N, 200, L, df, 20)
    end

    return nothing
end

main()