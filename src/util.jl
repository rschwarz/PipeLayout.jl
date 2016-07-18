# tolerance for numerical comparison
const TOL = 1e-6
const ɛ = TOL

nonzero(x::Real) = abs(x) > ɛ
