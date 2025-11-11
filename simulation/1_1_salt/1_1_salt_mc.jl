using Molly, PDMK4MC

n_atoms = 1000

charges = randn(n_atoms)
charges .-= sum(charges) / n_atoms

atoms = [Atom(mass=1.0u"g/mol", charge=charges[i], σ=0.5u"nm", ϵ=1.0u"kJ/mol") for i in 1:n_atoms]
boundary = CubicBoundary(30.0u"nm", 30.0u"nm", 30.0u"nm")

coords = place_atoms(n_atoms, boundary; min_dist=1.0u"nm")
pairwise_inters = (LennardJones(cutoff=DistanceCutoff(1.122462048309373u"nm"), shift = true),)

coords_val = zeros(3, n_atoms)
for i in 1:n_atoms
    for j in 1:3
        coords_val[j, i] = coords[i][j].val
    end
end

params = PDMK4MC.HPDMKParams(L = 30.0, digits = 3, n_per_leaf = 10, init = PDMK4MC.DIRECT)
tree = PDMK4MC.create_tree(coords_val, charges; params = params)

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

trial_args = Dict(:shift_size => 0.3u"nm")
sim = MetropolisMonteCarloPDMK(;
    temperature=300.0u"K",
    trial_moves=random_uniform_translation,
    trial_args=trial_args,
    tree=tree,
)

simulate!(sys, sim, 10_000)
PDMK4MC.destroy_tree!(tree)

println(sys.loggers.montecarlo.n_accept)

