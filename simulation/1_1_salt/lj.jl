using Molly
using Plots

inter_nonshifted = LennardJones(cutoff=DistanceCutoff(1.122462048309373u"nm"), shift = false)
inter_shifted = LennardJones(cutoff=DistanceCutoff(1.122462048309373u"nm"), shift = true)

sigma = 0.5u"nm"
epsilon = 2.0u"kJ/mol"

xs = range(0.5u"nm", 1.5u"nm", length=1000)
en = [Molly.pairwise_pe(inter_nonshifted, x, (sigma^2, epsilon)) for x in xs]
es = [Molly.pairwise_pe(inter_shifted, x, (sigma^2, epsilon)) for x in xs]

rxs = [x.val for x in xs]
ren = [e.val for e in en]
res = [e.val for e in es]

plot(rxs, ren, label="Nonshifted")
plot!(rxs, res, label="Shifted")
# ylims!(-0.05, 0.05)

plot(rxs, ren .- res, label="Nonshifted - Shifted")