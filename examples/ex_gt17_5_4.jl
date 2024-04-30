"""
	Example: ex_gt17_5_4.jl
	-------------------------------------------------------
	reference:
"""

include("../cim.jl")

function F(z::Complex)
    return [exp(im * z^2) 1.0; 1.0 1.0]
end

n = 200 # nbr of quadrature nodes
l = 1
r = 1
pbar = 6
circ = circle([0.0, 0.0], 3.0)
show_contr(circ)

pts = get_quadpts(circ, n)
show_quadpts(pts)

lambda = contr_int_ho(pts, F, 2, l, r, pbar)
show_eigs(lambda)