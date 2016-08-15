export ɛ

# tolerance for numerical comparison
const ɛ = 1e-6

nonzero(x::Real) = abs(x) > ɛ
