# extend Molly's MC process
using Molly
using Molly.Random
using PDMK4MC
using ParticleMeshEwald

export MetropolisMonteCarloPDMK, VerletFix
export random_uniform_translation, random_normal_translation, random_uniform_translation_fix

function random_uniform_translation(sys::Molly.System{D, <:Any, T}; shift_size=oneunit(eltype(eltype(sys.coords))), rng=Random.default_rng()) where {D, T}
    rand_idx = rand(rng, 1:length(sys.coords))
    direction = Molly.random_unit_vector(T, D, rng)
    magnitude = rand(rng, T) * shift_size
    return (rand_idx, magnitude * direction)
end 

# the first atom is fixed
function random_uniform_translation_fix(sys::Molly.System{D, <:Any, T}; shift_size=oneunit(eltype(eltype(sys.coords))), rng=Random.default_rng()) where {D, T}
    rand_idx = rand(rng, 2:length(sys.coords))
    direction = Molly.random_unit_vector(T, D, rng)
    magnitude = rand(rng, T) * shift_size
    return (rand_idx, magnitude * direction)
end 

function random_normal_translation(sys::Molly.System{D, <:Any, T}; shift_size=oneunit(eltype(eltype(sys.coords))), rng=Random.default_rng()) where {D, T}
    rand_idx = rand(rng, 1:length(sys.coords))
    direction = Molly.random_unit_vector(T, D, rng)
    magnitude = randn(rng, T) * shift_size
    return (rand_idx, magnitude * direction)
end

struct MetropolisMonteCarloPDMK{T, M, TE, TPME}
    temperature::T
    trial_moves::M
    trial_args::Dict
    eps_r::TE
    reconstruct::Int
    print_interval::Int
    energy_file
    check_accuracy::Bool
    n_check::Int
    pme::TPME
    accuracy_file
end

function MetropolisMonteCarloPDMK(; temperature, trial_moves, trial_args, eps_r, reconstruct, print_interval, energy_file, check_accuracy=false, n_check=100, pme=nothing, accuracy_file=nothing)
    return MetropolisMonteCarloPDMK(temperature, trial_moves, trial_args, eps_r, reconstruct, print_interval, energy_file, check_accuracy, n_check, pme, accuracy_file)
end

@inline function Molly.simulate!(sys::System,
                           sim::MetropolisMonteCarloPDMK,
                           n_steps::Integer, tree::PDMK4MC.Tree;
                           n_threads::Integer=Threads.nthreads(),
                           run_loggers=true,
                           rng=Random.default_rng())
    neighbors = Molly.find_neighbors(sys, sys.neighbor_finder; n_threads=n_threads)
    buffers = Molly.init_buffers!(sys, n_threads)
    E_old = Molly.potential_energy(sys, neighbors, buffers; n_threads=n_threads)
    PDMK4MC.form_outgoing_pw!(tree)
    PDMK4MC.form_incoming_pw!(tree)
    E_pdmk_old = PDMK4MC.eval_energy(tree) * 138.935457644u"kJ/mol" / sim.eps_r

    E_pme_old = 0.0
    E_pme_new = 0.0
    x = zeros(Float64, length(sys.coords))
    y = zeros(Float64, length(sys.coords))
    z = zeros(Float64, length(sys.coords))
    q = zeros(ComplexF64, length(sys.coords))

    for i in 1:length(sys.coords)
        q[i] = ComplexF64(sys.atoms[i].charge)
    end

    CSV.write(sim.accuracy_file, DataFrame(step_n = Int[], dE_pme = Float64[], dE_pdmk = Float64[], E_pme = Float64[]))
    CSV.write(sim.energy_file, DataFrame(step_n = Int[], E_lj = Float64[], E_elec = Float64[], E_total = Float64[]))

    for step_n in 1:n_steps

        if step_n % sim.reconstruct == 0
            tree_new = PDMK4MC.recontstruct_tree(tree)
            PDMK4MC.form_outgoing_pw!(tree_new)
            PDMK4MC.form_incoming_pw!(tree_new)
            E_pdmk_new = PDMK4MC.eval_energy(tree_new) * 138.935457644u"kJ/mol" / sim.eps_r

            with_logger(global_logger()) do
                @info "recontstructed, step_n = $(step_n), E_pdmk_old = $(E_pdmk_old), E_pdmk_new = $(E_pdmk_new)"
            end

            E_pdmk_old = E_pdmk_new
            PDMK4MC.destroy_tree!(tree)
            tree = tree_new
        end

        if sim.check_accuracy && (step_n % sim.n_check == 0)
            for i in 1:length(sys.coords)
                x[i] = sys.coords[i][1].val
                y[i] = sys.coords[i][2].val
                z[i] = sys.coords[i][3].val
            end
            E_pme_old = ParticleMeshEwald.energy(sim.pme, x, y, z, q)
        end

        idx, shift_vec = sim.trial_moves(sys; sim.trial_args...)
        # @show idx, shift_vec
        sys.coords[idx] = Molly.wrap_coords(sys.coords[idx] .+ shift_vec, sys.boundary)

        neighbors = Molly.find_neighbors(sys, sys.neighbor_finder; n_threads=n_threads)
        E_new = Molly.potential_energy(sys, neighbors, buffers, step_n; n_threads=n_threads)

        dE_pdmk = PDMK4MC.eval_shift_energy(tree, idx, shift_vec[1].val, shift_vec[2].val, shift_vec[3].val) * 138.935457644u"kJ/mol" / sim.eps_r

        if sim.check_accuracy && (step_n % sim.n_check == 0)
            x[idx] = sys.coords[idx][1].val
            y[idx] = sys.coords[idx][2].val
            z[idx] = sys.coords[idx][3].val
            E_pme_new = ParticleMeshEwald.energy(sim.pme, x, y, z, q)

            dE_pme = 4π * (E_pme_new - E_pme_old) * 138.935457644u"kJ/mol" / sim.eps_r
            CSV.write(sim.accuracy_file, DataFrame(step_n = step_n, dE_pme = dE_pme.val, dE_pdmk = dE_pdmk.val, E_pme = 4π * 138.935457644 * E_pme_new / sim.eps_r), append=true)
        end

        ΔE = E_new - E_old + dE_pdmk
        δ = ΔE / (sys.k * sim.temperature)
        # println("δ = $(δ), ΔE = $(ΔE), E_new = $(E_new), E_old = $(E_old), dE_pdmk = $(dE_pdmk)")
        if δ < 0 || (rand(rng) < exp(-δ))
            Molly.apply_loggers!(sys, nothing, neighbors, step_n, run_loggers; n_threads=n_threads,
                           current_potential_energy=(E_new + E_pdmk_old + dE_pdmk), success=true,
                           energy_rate=((E_new + E_pdmk_old + dE_pdmk) / (sys.k * sim.temperature)))
            PDMK4MC.update_shift!(tree, idx, shift_vec[1].val, shift_vec[2].val, shift_vec[3].val)
            E_old = E_new
            E_pdmk_old += dE_pdmk
        else
            sys.coords[idx] = Molly.wrap_coords(sys.coords[idx] .- shift_vec, sys.boundary)
            Molly.apply_loggers!(sys, nothing, neighbors, step_n, run_loggers; n_threads=n_threads,
                           current_potential_energy=E_old, success=false,
                           energy_rate=(E_old / (sys.k * sim.temperature)))
        end

        if step_n % sim.print_interval == 0
            CSV.write(sim.energy_file, DataFrame(step_n = step_n, E_lj = E_old.val, E_elec = E_pdmk_old.val, E_total = (E_old + E_pdmk_old).val), append=true)
        end

    end

    return sys
end

struct VerletFix{T, C}
    dt::T
    coupling::C
    remove_CM_motion::Int
    fix_idx::Int
end

function VerletFix(; dt, coupling=NoCoupling(), remove_CM_motion=1, fix_idx=1)
    return VerletFix(dt, coupling, Int(remove_CM_motion), fix_idx)
end

@inline function Molly.simulate!(sys,
                           sim::VerletFix,
                           n_steps::Integer;
                           n_threads::Integer=Threads.nthreads(),
                           run_loggers=true,
                           rng=Random.default_rng())
    needs_vir, needs_vir_steps = Molly.needs_virial_schedule(sim.coupling)
    sys.coords .= Molly.wrap_coords.(sys.coords, (sys.boundary,))
    !iszero(sim.remove_CM_motion) && Molly.remove_CM_motion!(sys)
    neighbors = Molly.find_neighbors(sys, sys.neighbor_finder; n_threads=n_threads)
    forces_t = Molly.zero_forces(sys)
    buffers = Molly.init_buffers!(sys, n_threads)
    Molly.apply_loggers!(sys, buffers, neighbors, 0, run_loggers; n_threads=n_threads)
    accels_t = forces_t ./ Molly.masses(sys)
    using_constraints = (length(sys.constraints) > 0)
    if using_constraints
        cons_coord_storage = zero(sys.coords)
    end

    for step_n in 1:n_steps
        needs_vir = (step_n % needs_vir_steps == 0)
        Molly.forces!(forces_t, sys, neighbors, buffers, Val(needs_vir), step_n; n_threads=n_threads)
        accels_t .= forces_t ./ Molly.masses(sys)

        sys.velocities .+= accels_t .* sim.dt

        sys.velocities[sim.fix_idx] = zero(sys.velocities[sim.fix_idx])

        if using_constraints
            cons_coord_storage .= sys.coords
        end
        sys.coords .+= sys.velocities .* sim.dt
        using_constraints && Molly.apply_position_constraints!(sys, cons_coord_storage;
                                                         n_threads=n_threads)

        if using_constraints
            sys.velocities .= (sys.coords .- cons_coord_storage) ./ sim.dt
        end

        sys.coords .= Molly.wrap_coords.(sys.coords, (sys.boundary,))

        if !iszero(sim.remove_CM_motion) && step_n % sim.remove_CM_motion == 0
            Molly.remove_CM_motion!(sys)
        end
        recompute_forces = Molly.apply_coupling!(sys, buffers, sim.coupling, sim, neighbors, step_n;
                                           n_threads=n_threads, rng=rng)

        neighbors = Molly.find_neighbors(sys, sys.neighbor_finder, neighbors, step_n, recompute_forces; n_threads=n_threads)

        Molly.apply_loggers!(sys, buffers, neighbors, step_n, run_loggers; n_threads=n_threads, current_forces=forces_t)
    end
    return sys
end