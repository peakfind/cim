"""
Example: prof_b12_4_11.jl
learn to use Profile and @profview in vscode
-------------------------------------------------------
reference: 
	
Example 4.11 in Wolf-Jurgen Beyn, An integral method for solving nonlinear eigenvalue problems
			
for use Profile in vscode, we refer to Julia in VS Code. What's New? | David Anthoff, Sebastian Pfitzner | JuliaCon 2022
https://www.youtube.com/watch?v=Okn_HKihWn8&t=0s
"""			

include("../cim.jl")
using Profile

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

# get the quadrature points
pts = get_quadpts(ellip, N)

lambda = @profview contr_int(pts, beyn12_411, 400, l)