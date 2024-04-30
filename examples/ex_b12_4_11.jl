"""
	Example: ex_b12_4_11.jl
	-------------------------------------------------------
	reference:
"""

include("../cim.jl")

ellip = ellipse([150.0, 0.0], 148.0, 148.0)
show_contr(ellip)

N = 50
pts = get_quadpts(ellip, N)
show_quadpts(pts)

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

lambda = contr_int(pts, beyn12_411, 400, 10)
show_eigs(lambda)