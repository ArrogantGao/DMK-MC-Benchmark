using PDMK4MC
using ParticleMeshEwald
using Random
using CSV, DataFrames
using BenchmarkTools

Random.seed!(2345)

# all particles are one a sphere of radius L / 4, centered as [L/2, L/2, L/2]

function benchmark_nonuniform_sys(N, ns, L, n_trials, df)

    coords = rand(Float64, 3, N)
    for i in 1:N
        xi = coords[1, i] - L / 2
        yi = coords[2, i] - L / 2
        zi = coords[3, i] - L / 2
        r = sqrt(xi^2 + yi^2 + zi^2)
        coords[1, i] = L / 2 + (L / 4) * xi / r
        coords[2, i] = L / 2 + (L / 4) * yi / r
        coords[3, i] = L / 2 + (L / 4) * zi / r
    end
    charges = [(-one(Float64))^i for i in 1:N]

    @info "Creating trees"
    tree_3 = PDMK4MC.create_tree(coords, charges; params=PDMK4MC.HPDMKParams(L = L, n_per_leaf=ns, init = DIRECT, digits = 3))
    @info "Tree 3 created"
    tree_6 = PDMK4MC.create_tree(coords, charges; params=PDMK4MC.HPDMKParams(L = L, n_per_leaf=ns, init = DIRECT, digits = 6))
    @info "Tree 6 created"
    tree_9 = PDMK4MC.create_tree(coords, charges; params=PDMK4MC.HPDMKParams(L = L, n_per_leaf=ns, init = DIRECT, digits = 9))
    @info "Tree 9 created"

    n_levels = tree_depth(tree_3)

    te_3 = 0.0
    tu_3 = 0.0
    te_6 = 0.0
    tu_6 = 0.0
    te_9 = 0.0
    tu_9 = 0.0

    for i in 1:n_trials

        idx = rand(1:N)
        dx = 0.25 * randn()
        dy = 0.25 * randn()
        dz = 0.25 * randn()

        te_3 += @elapsed PDMK4MC.eval_shift_energy(tree_3, idx, dx, dy, dz)
        tu_3 += @elapsed PDMK4MC.update_shift!(tree_3, idx, dx, dy, dz)
        te_6 += @elapsed PDMK4MC.eval_shift_energy(tree_6, idx, dx, dy, dz)
        tu_6 += @elapsed PDMK4MC.update_shift!(tree_6, idx, dx, dy, dz)
        te_9 += @elapsed PDMK4MC.eval_shift_energy(tree_9, idx, dx, dy, dz)
        tu_9 += @elapsed PDMK4MC.update_shift!(tree_9, idx, dx, dy, dz)
    end

    CSV.write(df, DataFrame(N = N, ns = ns, L = L, n_levels = n_levels, te_3 = te_3 / n_trials, tu_3 = tu_3 / n_trials, te_6 = te_6 / n_trials, tu_6 = tu_6 / n_trials, te_9 = te_9 / n_trials, tu_9 = tu_9 / n_trials), append=true)

    PDMK4MC.destroy_tree!(tree_3)
    PDMK4MC.destroy_tree!(tree_6)
    PDMK4MC.destroy_tree!(tree_9)
    GC.gc()
    return nothing
end

function main()
    df = joinpath(@__DIR__, "nonuniform_sys.csv")
    CSV.write(df, DataFrame(N = Int[], ns = Int[], L = Float64[], n_levels = Int[], te_3 = Float64[], tu_3 = Float64[], te_6 = Float64[], tu_6 = Float64[], te_9 = Float64[], tu_9 = Float64[]))
    
    rho = 1e4 / 64.63304^3

    N0 = 1000
    for i in 0:11
        N = N0 * 2^i
        L = (N / rho)^(1/3)
        benchmark_nonuniform_sys(N, 200, L, 100, df)
    end
    return nothing
end

main()