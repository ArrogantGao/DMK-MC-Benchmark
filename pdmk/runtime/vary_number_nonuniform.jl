using PDMK4MC
using CSV, DataFrames
using Random

Random.seed!(1234)


function main()

    filename = "pdmk_runtime_nonuniform_float.csv"
    CSV.write(filename, DataFrame(N = Int[], L = Float64[], ns = Int[], n_levels = Int[], time_propose = Float64[], time_accept = Float64[]))

    Ns = [10000, 12689, 16102, 20433, 25929, 32903, 41753, 52983, 67233, 85316, 108263, 137382, 174332, 221221, 280721, 356224, 452036, 573615, 727895, 923670, 1172100, 1487352, 1887392, 2395026, 3039196, 3856622, 4893900, 6210168, 7880462, 10000000]
    nss = [200]

    for N in Ns
        for ns in nss
            println("Running for N = $N, ns = $ns")
            L = Float32(3.0 * ( N^(1.0/3.0)))

            coords = rand(Float32, 3, N)
            for i in 1:N
                xi = coords[1, i] - L / 2
                yi = coords[2, i] - L / 2
                zi = coords[3, i] - L / 2
                r = sqrt(xi^2 + yi^2 + zi^2)
                coords[1, i] = L / 2 + (L / 4) * xi / r
                coords[2, i] = L / 2 + (L / 4) * yi / r
                coords[3, i] = L / 2 + (L / 4) * zi / r
            end
            charges = [(-one(Float32))^i for i in 1:N]

            params = PDMK4MC.HPDMKParams(L = L, n_per_leaf=ns, init = PROXY)
            tree = PDMK4MC.create_tree(coords, charges; params=params)

            println("Tree created")

            n_levels = tree_depth(tree)

            # warm up
            eval_shift_energy(tree, 1, 0.1f0, 0.1f0, 0.1f0)
            update_shift!(tree, 1, 0.1f0, 0.1f0, 0.1f0)

            println("Warming up done")

            time_propose = 0.0
            time_accept = 0.0

            n_trials = 10000
            for i in 1:n_trials

                idx = rand(1:N)

                tx_new = rand(Float32) - 0.5f0
                ty_new = rand(Float32) - 0.5f0
                tz_new = rand(Float32) - 0.5f0
                r_new = sqrt(tx_new^2 + ty_new^2 + tz_new^2)

                x_new = L / 2 + (L / 4) * tx_new / r_new
                y_new = L / 2 + (L / 4) * ty_new / r_new
                z_new = L / 2 + (L / 4) * tz_new / r_new

                dx = x_new - coords[1, idx]
                dy = y_new - coords[2, idx]
                dz = z_new - coords[3, idx]

                time_propose += @elapsed eval_shift_energy(tree, idx, dx, dy, dz)
                time_accept += @elapsed update_shift!(tree, idx, dx, dy, dz)
            end

            destroy_tree!(tree)

            println("Time propose: $(time_propose / n_trials), Time accept: $(time_accept / n_trials)")

            CSV.write(filename, DataFrame(N = N, L = L, ns = ns, n_levels = n_levels,time_propose = time_propose / n_trials, time_accept = time_accept / n_trials), append=true)
        end
    end
end

main()