"""
	Example: ex_neppack_tutorial2.jl
	-------------------------------------------------------
	reference: Tutorials for NEP-PACK Tutorial 2 (Contour)
"""

include("../cim.jl")
using BenchmarkTools

# define the matrix of the NEP
function NEP(z::Complex)
	n = 1000 
	A0 = diagm(0 => ones(n))
	A1 = diagm(-2 => ones(n-2), 0 => 30 * (n:-1:1)/n, 1 => 3*ones(n-1))/3
	A2 = diagm(-1 => ones(n-1), 0 => (1:n)/n, 1 => sin.(range(0,5,length=n-1)))/10

	return A0 + A1 * z + A2 * exp(-z)
end

# parameters
N = 25 # the number of the quadrature nodes
l = 10
r = 10
pbar = 6

circ = circle([0.0, 0.0], 0.5)
pts = get_quadpts(circ, N)

# lambda = @btime contr_int(pts, NEP, 400, l)
lambda = contr_int_ho(pts, NEP, 1000, l, r, pbar)
# @show lambda