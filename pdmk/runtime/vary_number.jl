using PDMK4MC
using CSV, DataFrames
using Random

Random.seed!(1234)


function main()

    filename = "pdmk_runtime_float.csv"
    CSV.write(filename, DataFrame(N = Int[], L = Float64[], time_propose = Float64[], time_accept = Float64[]))

    Ns = [10000, 12689, 16102, 20433, 25929, 32903, 41753, 52983, 67233, 85316, 108263, 137382, 174332, 221221, 280721, 356224, 452036, 573615, 727895, 923670, 1172102, 1487352, 1887391, 2395026, 3039195, 3856620, 4893900, 6210169, 7880462, 10000000]

    for N in Ns
        println("Running for N = $N")
        L = Float32(3.0 * ( N^(1.0/3.0)))

        coords = rand(Float32, 3, N) .* L
        charges = randn(Float32, N)
        charges .-= sum(charges) / N

        params = PDMK4MC.HPDMKParams(L = L, n_per_leaf=200, init = PROXY)
        tree = PDMK4MC.create_tree(coords, charges; params=params)

        # warm up
        eval_shift_energy(tree, 1, 0.1f0, 0.1f0, 0.1f0)
        update_shift!(tree, 1, 0.1f0, 0.1f0, 0.1f0)

        time_propose = 0.0
        time_accept = 0.0

        n_trials = 10000
        for i in 1:n_trials

            idx = rand(1:N)
            dx = Float32(0.25 * randn())
            dy = Float32(0.25 * randn())
            dz = Float32(0.25 * randn())

            time_propose += @elapsed eval_shift_energy(tree, idx, dx, dy, dz)
            time_accept += @elapsed update_shift!(tree, idx, dx, dy, dz)
        end

        destroy_tree!(tree)

        println("Time propose: $(time_propose / n_trials), Time accept: $(time_accept / n_trials)")

        CSV.write(filename, DataFrame(N = N, L = L, time_propose = time_propose / n_trials, time_accept = time_accept / n_trials), append=true)
    end
end

main()