using Molly
using LinearAlgebra

using Molly: force_cutoff, pe_cutoff

export ShiftedLennardJones

@kwdef struct ShiftedLennardJones{C, H, S, E, W} <: Molly.PairwiseInteraction
    cutoff::C = NoCutoff()
    use_neighbors::Bool = false
    shortcut::H = Molly.lj_zero_shortcut
    σ_mixing::S = Molly.lorentz_σ_mixing
    ϵ_mixing::E = Molly.geometric_ϵ_mixing
    weight_special::W = 1
    shift::Bool = false
end

use_neighbors(inter::ShiftedLennardJones) = inter.use_neighbors

function Base.zero(lj::ShiftedLennardJones{C, H, S, E, W}) where {C, H, S, E, W}
    return ShiftedLennardJones(
        lj.cutoff,
        lj.use_neighbors,
        lj.shortcut,
        lj.σ_mixing,
        lj.ϵ_mixing,
        zero(W),
    )
end

function Base.:+(l1::ShiftedLennardJones, l2::ShiftedLennardJones)
    return ShiftedLennardJones(
        l1.cutoff,
        l1.use_neighbors,
        l1.shortcut,
        l1.σ_mixing,
        l1.ϵ_mixing,
        l1.weight_special + l2.weight_special,
    )
end

function Molly.inject_interaction(inter::ShiftedLennardJones, params_dic)
    key_prefix = "inter_LJ_"
    return ShiftedLennardJones(
        inter.cutoff,
        inter.use_neighbors,
        inter.shortcut,
        inter.σ_mixing,
        inter.ϵ_mixing,
        dict_get(params_dic, key_prefix * "weight_14", inter.weight_special),
    )
end

function Molly.extract_parameters!(params_dic, inter::ShiftedLennardJones, ff)
    key_prefix = "inter_LJ_"
    params_dic[key_prefix * "weight_14"] = inter.weight_special
    return params_dic
end

@inline function Molly.force(inter::ShiftedLennardJones,
                       dr,
                       atom_i,
                       atom_j,
                       force_units=u"kJ * mol^-1 * nm^-1",
                       special=false,
                       args...)
    if inter.shortcut(atom_i, atom_j)
        return ustrip.(zero(dr)) * force_units
    end
    σ = inter.σ_mixing(atom_i, atom_j)
    ϵ = inter.ϵ_mixing(atom_i, atom_j)

    cutoff = inter.cutoff
    r = norm(dr)
    σ2 = σ^2
    params = (σ2, ϵ)

    f = force_cutoff(cutoff, inter, r, params)
    fdr = (f / r) * dr
    if special
        return fdr * inter.weight_special
    else
        return fdr
    end
end

function Molly.pairwise_force(::ShiftedLennardJones, r, (σ2, ϵ))
    six_term = (σ2 / r^2) ^ 3
    return (24ϵ / r) * (2 * six_term ^ 2 - six_term)
end

@inline function Molly.potential_energy(inter::ShiftedLennardJones,
                                  dr,
                                  atom_i,
                                  atom_j,
                                  energy_units=u"kJ * mol^-1",
                                  special=false,
                                  args...)
    if inter.shortcut(atom_i, atom_j)
        return ustrip(zero(dr[1])) * energy_units
    end
    σ = inter.σ_mixing(atom_i, atom_j)
    ϵ = inter.ϵ_mixing(atom_i, atom_j)

    cutoff = inter.cutoff
    r = norm(dr)
    σ2 = σ^2
    params = (σ2, ϵ)

    pe = Molly.pe_cutoff(cutoff, inter, r, params)
    if special
        return pe * inter.weight_special
    else
        return pe
    end
end

function Molly.pairwise_pe(inter::ShiftedLennardJones, r, (σ2, ϵ))
    six_term = (σ2 / r^2) ^ 3
    if inter.shift
        six_term_cutoff = (σ2 / inter.cutoff.dist_cutoff^2) ^ 3
        return (4ϵ * (six_term ^ 2 - six_term) - 4ϵ * (six_term_cutoff ^ 2 - six_term_cutoff)) * (r <= inter.cutoff.dist_cutoff)
    else
        return 4ϵ * (six_term ^ 2 - six_term)
    end
end