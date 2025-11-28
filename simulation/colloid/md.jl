using Molly, PDMK4MC
using JLD2, CSV, DataFrames
using Random
using LoggingExtras

Random.seed!(1234)
logger = SimpleLogger(open(joinpath(@__DIR__, "md.log"), "w+"))
global_logger(logger)

include(joinpath(@__DIR__, "../MollyExt/simulator.jl"))
include(joinpath(@__DIR__, "../MollyExt/truncated_lj.jl"))

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

    pairwise_inters = (ShiftedLennardJones(cutoff=DistanceCutoff(3.0u"nm"), use_neighbors = true, shift = true), CoulombEwald(dist_cutoff=3.0u"nm", coulomb_const = Molly.coulomb_const / eps_r, use_neighbors=true))
    general_inters = (Molly.PME(3.0u"nm", atoms, boundary, ϵr = eps_r),)

    sys = System(
        atoms=atoms,
        coords=coords,
        boundary=boundary,
        pairwise_inters=pairwise_inters,
        general_inters=general_inters,
        loggers=(
            coords=CoordinatesLogger(100, dims=3),
            temp=TemperatureLogger(100),
            energy=PotentialEnergyLogger(100),
        ),
        neighbor_finder=DistanceNeighborFinder(eligible=trues(n_atoms, n_atoms), dist_cutoff=3.0u"nm", n_steps=100),
    )

    minimizer = SteepestDescentMinimizer()
    simulate!(sys, minimizer)
    coords[1] = Molly.SVector(L_box / 2.0, L_box / 2.0, L_box / 2.0)

    random_velocities!(sys, temp)
    sys.velocities[1] = zero(sys.velocities[1])

    simulator = VerletFix(dt=0.0001u"ps", coupling = AndersenThermostat(temp, 1.0u"ps"), fix_idx=1)
    
    simulate!(sys, simulator, 100_000)

    simulator = VerletFix(dt=0.001u"ps", coupling = AndersenThermostat(temp, 1.0u"ps"), fix_idx=1)
    simulate!(sys, simulator, 1_000_000)

    jldsave(joinpath(@__DIR__, "data/md.jld2"), sys=sys)
end

main()