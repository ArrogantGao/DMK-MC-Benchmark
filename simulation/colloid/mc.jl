using Molly, PDMK4MC
using JLD2, CSV, DataFrames
using Random
using LoggingExtras

Random.seed!(1234)

include(joinpath(@__DIR__, "../MollyExt/simulator.jl"))
include(joinpath(@__DIR__, "../MollyExt/truncated_lj.jl"))
logger = SimpleLogger(open(joinpath(@__DIR__, "mc.log"), "w+"))
global_logger(logger)

function main()
    L_box = 9.4u"nm" # 1mol/L NaCl in water

    eps_r = 78.4 # water

    σ_Na = 0.1369u"nm"
    eps_Na = 0.0874393 * 4.184u"kJ/mol"
    σ_Cl = 0.2513u"nm"
    eps_Cl = 0.0355910 * 4.184u"kJ/mol"

    temp = 300.0u"K"

    atom_colloid = [Atom(mass=(1.0)u"g/mol", charge=+100.0, σ=1.0u"nm", ϵ=eps_Na)]
    atoms_na = [Atom(mass=1.0u"g/mol", charge=+1.0, σ=σ_Na, ϵ=eps_Na) for _ in 1:400]
    atoms_cl = [Atom(mass=1.0u"g/mol", charge=-1.0, σ=σ_Cl, ϵ=eps_Cl) for _ in 1:500]
    atoms = vcat(atom_colloid, atoms_na, atoms_cl)
    boundary = CubicBoundary(L_box, L_box, L_box)

    n_atoms = length(atoms)

    coords = place_atoms(n_atoms, boundary; min_dist=0.3u"nm")
    coords[1] = Molly.SVector(L_box / 2.0, L_box / 2.0, L_box / 2.0)

    pairwise_inters = (ShiftedLennardJones(cutoff=DistanceCutoff(3.0u"nm"), shift = true, use_neighbors = true),)

    coords_val = zeros(3, n_atoms)
    for i in 1:n_atoms
        for j in 1:3
            coords_val[j, i] = coords[i][j].val
        end
    end
    charges_val = [atom.charge for atom in atoms]

    params = PDMK4MC.HPDMKParams(L = L_box.val, digits = 3, n_per_leaf = 50, init = PDMK4MC.DIRECT)
    tree = PDMK4MC.create_tree(coords_val, charges_val; params = params)

    sys = System(
        atoms=atoms,
        coords=coords,
        boundary=boundary,
        pairwise_inters=pairwise_inters,
        neighbor_finder=DistanceNeighborFinder(eligible=trues(n_atoms, n_atoms), dist_cutoff=3.0u"nm", n_steps=100),
        loggers=(
            coords=CoordinatesLogger(1000, dims=3),
            montecarlo=MonteCarloLogger(),
        ),
    )

    trial_args = Dict(:shift_size => 0.3u"nm")
    sim = MetropolisMonteCarloPDMK(;
        temperature=300.0u"K",
        trial_moves=random_uniform_translation_fix,
        trial_args=trial_args,
        eps_r=eps_r,
        reconstruct = 5000,
        print_interval = 100,
        energy_file = joinpath(@__DIR__, "data/energy_mc_accuracy.csv"),
        check_accuracy = true,
        n_check = 1000,
        pme = ParticleMeshEwald.PME(2.0, (9.4, 9.4, 9.4), 3.0, n_atoms),
        accuracy_file = joinpath(@__DIR__, "data/dE_mc_accuracy.csv")
    )

    simulate!(sys, sim, 10_000_000, tree)
    jldsave(joinpath(@__DIR__, "data/mc.jld2"), sys=sys)
end

main()