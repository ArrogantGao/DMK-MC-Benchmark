using Molly, PDMK4MC
using Random

Random.seed!(1234)

include(joinpath(@__DIR__, "../molly.jl"))
# logger = SimpleLogger(open(joinpath(@__DIR__, "NaCl_mc.log"), "w+"))
# global_logger(logger)

function main()
    n_atoms = 1000
    L_box = 9.4u"nm" # 1mol/L NaCl in water

    eps_r = 78.4 # water

    σ_Na = 0.1369u"nm"
    eps_Na = 0.0874393 * 4.184u"kJ/mol"
    σ_Cl = 0.2513u"nm"
    eps_Cl = 0.0355910 * 4.184u"kJ/mol"

    atoms_na = [Atom(mass=1.0u"g/mol", charge=+1.0, σ=σ_Na, ϵ=eps_Na) for _ in 1:n_atoms / 2]
    atoms_cl = [Atom(mass=1.0u"g/mol", charge=-1.0, σ=σ_Cl, ϵ=eps_Cl) for _ in 1:n_atoms / 2]
    atoms = vcat(atoms_na, atoms_cl)
    boundary = CubicBoundary(L_box, L_box, L_box)

    coords = place_atoms(n_atoms, boundary; min_dist=0.2u"nm")
    pairwise_inters = (LennardJones(cutoff=DistanceCutoff(3.0u"nm"), shift = true),)

    coords_val = zeros(3, n_atoms)
    for i in 1:n_atoms
        for j in 1:3
            coords_val[j, i] = coords[i][j].val
        end
    end
    charges_val = [atom.charge for atom in atoms]

    params = PDMK4MC.HPDMKParams(L = L_box.val, digits = 6, n_per_leaf = 30, init = PDMK4MC.DIRECT)
    tree = PDMK4MC.create_tree(coords_val, charges_val; params = params)

    energy_dmk = PDMK4MC.eval_energy(tree) * 138.935457644u"kJ/mol" / eps_r

    sys = System(
        atoms=atoms,
        coords=coords,
        boundary=boundary,
        pairwise_inters=pairwise_inters,
        loggers=(
            coords=CoordinatesLogger(n_atoms, dims=3),
            montecarlo=MonteCarloLogger(),
        ),
    )

    trial_args = Dict(:shift_size => 0.1u"nm")
    sim = MetropolisMonteCarloPDMK(;
        temperature=300.0u"K",
        trial_moves=random_uniform_translation,
        trial_args=trial_args,
        eps_r=eps_r,
        reconstruct = 5000,
        print_interval = 100
    )

    simulate!(sys, sim, 500_000, tree)
end