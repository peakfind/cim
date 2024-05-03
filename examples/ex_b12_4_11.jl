"""
	Example: ex_b12_4_11.jl
	-------------------------------------------------------
	reference: Example 4.11 in Wolf-Jurgen Beyn, An integral 
	           method for solving nonlinear eigenvalue problems
"""

include("../cim.jl")
using BenchmarkTools

# define the matrix of the NEP
function beyn12_411(z::Complex)
    D = 400
    NEP = zeros(ComplexF64, D, D)
    ed = 2.0*D - 4.0*z/(6.0*D) # diagonal entry
	eoffd = -1.0*D - z/(6.0*D) # off-diagonal entry

	# top row
	NEP[1, 1] = ed
	NEP[1, 2] = eoffd

	# interior rows
	for d = 2:D-1
		NEP[d,d] = ed
		NEP[d,d-1] = NEP[d,d+1] = eoffd
    end

	NEP[D, D-1] = eoffd # bottom row
	NEP[D, D] = 0.5*ed + z/(z-1.0)

	return NEP
end

# parameters
N = 50 # the number of the quadrature nodes
l = 10

# define the contour
ellip = ellipse([150.0, 0.0], 148.0, 148.0)
show_contr(ellip)

# get the quadrature points
pts = get_quadpts(ellip, N)
show_quadpts(pts)

lambda = contr_int(pts, beyn12_411, 400, l)
show_eigs(lambda)