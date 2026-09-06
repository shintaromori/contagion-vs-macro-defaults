module Quadrature

using FastGaussQuadrature

export get_hermite_quadrature

"""
    get_hermite_quadrature(GH_N::Int)

Returns Gauss-Hermite nodes `x` and normalized log-weights `log_w` for `GH_N` quadrature points.
Matches `numpy.polynomial.hermite.hermgauss(GH_N)` convention:
- `w_nodes = w / sqrt(pi)`
- `w_nodes = w_nodes / sum(w_nodes)`
- `log_w = log.(w_nodes)`
"""
function get_hermite_quadrature(GH_N::Int)
    x, w = gausshermite(GH_N)
    w_norm = w ./ sqrt(pi)
    w_norm ./= sum(w_norm)
    log_w = log.(w_norm)
    return x, log_w
end

end # module Quadrature
