using LinearAlgebra

function charge_integration(coords, charges, center, bins, L)
    N = length(charges)
    cbin = zeros(bins)
    ci = zeros(bins)
    d = L / bins
    for coord in coords
        for i in 1:N
            for mx in -1:1, my in -1:1, mz in -1:1
                rs = norm(coord[i] .- center .+ L * [mx, my, mz])
                if rs < L
                    idx = Int(floor(rs / d)) + 1
                    cbin[idx] += charges[i]
                end
            end
        end
    end
    cbin ./= length(coords)

    ci[1] = cbin[1]
    for i in 2:bins
        ci[i] = ci[i - 1] + cbin[i]
    end

    return ci
end