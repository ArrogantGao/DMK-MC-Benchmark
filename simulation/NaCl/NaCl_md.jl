using Molly, PDMK4MC
using JLD2, CSV, DataFrames
using Random
using LoggingExtras

Random.seed!(1234)
logger = SimpleLogger(open(joinpath(@__DIR__, "NaCl_md.log"), "w+"))
global_logger(logger)


function main()
    n_atoms = 1000
    L_box = 9.4u"nm" # 1mol/L NaCl in water

    eps_r = 78.4 # water

    σ_Na = 0.1369u"nm"
    eps_Na = 0.0874393 * 4.184u"kJ/mol"
    σ_Cl = 0.2513u"nm"
    eps_Cl = 0.0355910 * 4.184u"kJ/mol"

    temp = 300.0u"K"

    atoms_na = [Atom(mass=1.0u"g/mol", charge=+1.0, σ=σ_Na, ϵ=eps_Na) for _ in 1:n_atoms / 2]
    atoms_cl = [Atom(mass=1.0u"g/mol", charge=-1.0, σ=σ_Cl, ϵ=eps_Cl) for _ in 1:n_atoms / 2]
    atoms = vcat(atoms_na, atoms_cl)
    boundary = CubicBoundary(L_box, L_box, L_box)

    coords = place_atoms(n_atoms, boundary; min_dist=0.2u"nm")

    pairwise_inters = (LennardJones(cutoff=DistanceCutoff(3.0u"nm"), use_neighbors = true), CoulombEwald(dist_cutoff=3.0u"nm", coulomb_const = Molly.coulomb_const / eps_r, use_neighbors=true))
    general_inters = (PME(3.0u"nm", atoms, boundary, ϵr = eps_r),)

    sys = System(
        atoms=atoms,
        coords=coords,
        boundary=boundary,
        pairwise_inters=pairwise_inters,
        general_inters=general_inters,
        loggers=(
            coords=CoordinatesLogger(n_atoms, dims=3),
            temp=TemperatureLogger(100),
            energy=PotentialEnergyLogger(100),
        ),
        neighbor_finder=DistanceNeighborFinder(eligible=trues(n_atoms, n_atoms), dist_cutoff=3.0u"nm", n_steps=100),
    )

    random_velocities!(sys, temp)

    simulator = Verlet(dt=0.001u"ps", coupling = AndersenThermostat(temp, 1.0u"ps"))
    
    simulate!(sys, simulator, 500_000)

    jldsave(joinpath(@__DIR__, "NaCl_md.jld2"), sys=sys)
end

main()