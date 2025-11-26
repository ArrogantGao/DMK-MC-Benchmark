using Molly

function rdf_type(coords_1, coords_2, L, npoints=100)
    n_1 = length(coords_1)
    n_2 = length(coords_2)

    rs = typeof(L)[]
    for i in 1:n_1
        for j in 1:n_2
            for mx in -1:1, my in -1:1, mz in -1:1
                drx = coords_1[i][1] - coords_2[j][1] + L * mx
                dry = coords_1[i][2] - coords_2[j][2] + L * my
                drz = coords_1[i][3] - coords_2[j][3] + L * mz
                dr = sqrt(drx^2 + dry^2 + drz^2)
                if dr > zero(L) && dr < L
                    push!(rs, dr)
                end
            end
        end
    end

    d = L / npoints
    bins = range(zero(L), L, length=npoints + 1)
    counts = zeros(Int, npoints)
    for r in rs
        idx = Int(floor(r / d)) + 1
        @assert bins[idx] <= r < bins[idx + 1]
        counts[idx] += 1
    end

    density = [counts[i] / n_1 / (4π * ((bins[i] + bins[i + 1]) / 2)^2 * d) for i in 1:npoints]

    return bins, density
end