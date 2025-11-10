using Molly

n_atoms = 330
atoms = [Atom(mass=1.0u"g/mol", charge=0.0, σ=0.5u"nm", ϵ=1.0u"kJ/mol") for i in 1:n_atoms]
boundary = CubicBoundary(30.0u"nm", 30.0u"nm", 30.0u"nm")

coords = place_atoms(n_atoms, boundary; min_dist=1.0u"nm")
pairwise_inters = (LennardJones(cutoff=DistanceCutoff(1.122462048309373u"nm"), shift = true),)

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
sim = MetropolisMonteCarlo(;
    temperature=300.0u"K",
    trial_moves=random_uniform_translation!,
    trial_args=trial_args,
)

simulate!(sys, sim, 10_000)

println(sys.loggers.montecarlo.n_accept)

