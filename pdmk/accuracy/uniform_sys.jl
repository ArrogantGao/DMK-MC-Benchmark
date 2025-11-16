using PDMK4MC
using ParticleMeshEwald
using Random
using CSV, DataFrames

Random.seed!(2345)

function accuracy_uniform_sys(N, ns, L, n_trials, df)

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

    # alpha = 0.35 # magic number to balance short and long
    # pme = PME(alpha, (L, L, L), 5.0, N)
    # @info "PME created"
    # E_old = ParticleMeshEwald.energy(pme, x, y, z, q)
    # @info "PME Energy evaluated"

    for i in 1:n_trials

        idx = rand(1:N)
        dx = 0.25 * randn()
        dy = 0.25 * randn()
        dz = 0.25 * randn()

        dE_dmk_3 = PDMK4MC.eval_shift_energy(tree_3, idx, dx, dy, dz) / 4π
        dE_dmk_6 = PDMK4MC.eval_shift_energy(tree_6, idx, dx, dy, dz) / 4π
        dE_dmk_9 = PDMK4MC.eval_shift_energy(tree_9, idx, dx, dy, dz) / 4π
        dE_dmk_12 = PDMK4MC.eval_shift_energy(tree_12, idx, dx, dy, dz) / 4π

        # x_old, y_old, z_old = x[idx], y[idx], z[idx]
        # x[idx] = (x_old + dx + L) % L
        # y[idx] = (y_old + dy + L) % L
        # z[idx] = (z_old + dz + L) % L
        # E_new = ParticleMeshEwald.energy(pme, x, y, z, q)

        abs_err_3 = abs(dE_dmk_3 - dE_dmk_12)
        abs_err_6 = abs(dE_dmk_6 - dE_dmk_12)
        abs_err_9 = abs(dE_dmk_9 - dE_dmk_12)
        rel_err_3 = abs(dE_dmk_3 - dE_dmk_12) / abs(dE_dmk_12)
        rel_err_6 = abs(dE_dmk_6 - dE_dmk_12) / abs(dE_dmk_12)
        rel_err_9 = abs(dE_dmk_9 - dE_dmk_12) / abs(dE_dmk_12)

        # x[idx] = x_old
        # y[idx] = y_old
        # z[idx] = z_old

        # CSV.write(df, DataFrame(N = N, ns = ns, L = L, digits = digits, n_trials = n_trials, n_levels = n_levels, i = i, abs_err_3 = abs_err_3, rel_err_3 = rel_err_3, abs_err_6 = abs_err_6, rel_err_6 = rel_err_6, abs_err_9 = abs_err_9, rel_err_9 = rel_err_9, abs_err_12 = abs_err_12, rel_err_12 = rel_err_12), append=true)
        CSV.write(df, DataFrame(N = N, ns = ns, L = L, n_trials = n_trials, n_levels = n_levels, i = i, abs_err_3 = abs_err_3, rel_err_3 = rel_err_3, abs_err_6 = abs_err_6, rel_err_6 = rel_err_6, abs_err_9 = abs_err_9, rel_err_9 = rel_err_9), append=true)
    end

    PDMK4MC.destroy_tree!(tree_3)
    PDMK4MC.destroy_tree!(tree_6)
    PDMK4MC.destroy_tree!(tree_9)
    PDMK4MC.destroy_tree!(tree_12)
    GC.gc()
    return nothing
end

function main()
    df = joinpath(@__DIR__, "uniform_sys.csv")
    CSV.write(df, DataFrame(N = Int[], ns = Int[], L = Float64[], n_trials = Int[], n_levels = Int[], i = Int[], abs_err_3 = Float64[], rel_err_3 = Float64[], abs_err_6 = Float64[], rel_err_6 = Float64[], abs_err_9 = Float64[], rel_err_9 = Float64[]))
    
    rho = 1e4 / 64.63304^3

    accuracy_uniform_sys(1000, 50, (1000/rho)^(1/3), 20, df)

    N0 = 2000
    for i in 0:10
        N = N0 * 2^i
        L = (N / rho)^(1/3)
        accuracy_uniform_sys(N, 200, L, 20, df)
    end
    return nothing
end

main()