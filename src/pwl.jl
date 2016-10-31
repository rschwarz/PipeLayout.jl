# utilities for piece-wise linear functions

"slope and offset for line through two points"
function linecoefs(x1, y1, x2, y2)
    x1 ≈ x2 && error("line too steep!")

    # line equation:             y = slope*x + offset
    # from Taylor of 1st order:      f'(x0)*(x - x0) + f(x0)
    # use x1 for x0:             y = slope*(x - x1) + f(x1)
    slope = (y2 - y1)/(x2 - x1)
    offset = y1 - slope*x1
    slope, offset
end

"""
Return coefficient matrix [a b c] of overestimating inequalities of the form
  a*x + b*y ≥ c
for a piece-wise linear function given by
"""
function pwl_ineqs(xs, ys)
    @assert length(xs) == length(ys)
    m = length(xs) - 1

    ineqs = fill(0.0, (m, 3))
    for i = 1:m
        slope, offset = linecoefs(xs[i], ys[i], xs[i+1], ys[i+1])
        ineqs[i, :] = [-slope, 1.0, offset]
    end
    ineqs
end

"Compute the inverse of a monotone PWL function at a point"
function pwl_inverse(xs, ys, y)
    @assert length(xs) == length(ys)
    @assert issorted(xs)
    increasing, decreasing = issorted(ys), issorted(ys, rev=true)
    @assert increasing || decreasing
    if decreasing # make it increasing
        y, ys = -y, -ys
    end
    @assert ys[1] ≤ y ≤ ys[end]

    right = findfirst(yi -> y ≤ yi, ys)
    if right == 1 # on the left boundary
        return xs[1]
    end

    left = right - 1
    @assert ys[left] ≤ y ≤ ys[right]
    if ys[left] ≈ ys[right]
        return xs[left]
    end
    # y = (1-λ)*ys[left] + λ*ys[right]
    λ = (y - ys[left])/(ys[right] - ys[left])
    x = (1-λ)*xs[left] + λ*xs[right]
    @assert xs[left] ≤ x ≤ xs[right]
    x
end
