# utilities for pipe diameters

"compute equivalent diameter from serial pipe segments"
function serialmerge(diams::Vector{Diameter}, lengths::Vector{Float64})
    @assert length(diams) == length(lengths)
    @assert sum(lengths) ≈ 1.0
    y = sum(l / d.value^5 for (d,l) in zip(diams, lengths))
    y^(-1/5)
end

"compute lengths of two segments for equivalent diameter"
function pipesplit(diams::Vector, equiv::Float64)
    ds = [d.value for d in diams]
    dmin, dmax = extrema(ds)
    @assert dmin ≈ equiv || dmin < equiv < dmax || equiv ≈ dmax
    y, ys = equiv^(-5), ds.^(-5)
    @assert issorted(ys, rev=true) # decreasing

    l = fill(0.0, length(diams))

    right = findfirst(yi -> y ≥ yi, ys) # decreasing
    if right == 1 # on left boundary
        l[1] = 1.0
        return l
    end

    left = right - 1
    λ = (y - ys[left])/(ys[right] - ys[left])
    l[left] = 1 - λ
    l[right] = λ
    l
end
