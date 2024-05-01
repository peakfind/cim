"""
	Example: ex_gt17_5_4.jl
	-------------------------------------------------------
	reference: Section 5.4 Numerical illustration
	           Stefan Guttel, Francoise Tisseur, The nonlinear 
			   eigenvalue problem, Acta Numer., 2017.
"""

include("../cim.jl")

# define the matrix of the NEP
function F(z::Complex)
    return [exp(im * z^2) 1.0; 1.0 1.0]
end

# parameters
n    = 200 # number of the quadrature nodes
l    = 1
r    = 1
pbar = 6

# define the contour
circ = circle([0.0, 0.0], 3.0)
show_contr(circ)

# get the quadrature points
pts = get_quadpts(circ, n)
show_quadpts(pts)

lambda = contr_int_ho(pts, F, 2, l, r, pbar)
show_eigs(lambda)