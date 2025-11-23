using Molly, PDMK4MC

n_atoms = 500

charges = randn(n_atoms)
charges .-= sum(charges) / n_atoms

sum(charges)

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

params = PDMK4MC.HPDMKParams(L = 30.0, digits = 6, n_per_leaf = 10, init = PDMK4MC.DIRECT)
tree = PDMK4MC.create_tree(coords_val, charges; params = params)

energy_dmk = PDMK4MC.eval_energy(tree) * 138.935457644u"kJ/mol"


d = 5.0u"nm"
ewald_short = CoulombEwald(dist_cutoff = d, use_neighbors = true)
ewald_long = Ewald(d)

sys = System(
    atoms=atoms,
    coords=coords,
    boundary=boundary,
    # pairwise_inters=(ewald_short,),
    # general_inters=(ewald_long,),
    pairwise_inters=(LennardJones(cutoff=DistanceCutoff(d), shift = true),),
    neighbor_finder=DistanceNeighborFinder(
        eligible=trues(n_atoms, n_atoms),
        n_steps=10,
        dist_cutoff=d,
    ),
)

energy_ewald = potential_energy(sys)


PDMK4MC.destroy_tree!(tree)