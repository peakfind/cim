#=---------------------------------------------------------
- cim.jl
- TODO: we should write a script to implement the method first
-
---------------------------------------------------------=#

using LinearAlgebra
using Plots

# contours that we plan to support
abstract type AbstractContour end

# TODO: Can we provide a generic method for different types of contours?
struct ellipse <: AbstractContour
    center::ComplexF64 # TODO: store the vector{F64} or complex?
    semi_x::Float64
    semi_y::Float64
end

struct rectangle <: AbstractContour
    center::Vector{Float64}
    semi_x::Float64
    semi_y::Float64
end

struct circle <: AbstractContour
    center::Vector{Float64}
    radius::Float64
end

function show_contr(contour::AbstractContour)
    if typeof(contour) == ellipse
        # plot()
    elseif typeof(contour) == circle
        # plot()
    elseif typeof(contour) == rectangle
        # plot()
    else
        # throw some error : unsupported type of contour
        throw(ArgumentError("unspported type of contour"))
    end
end

# TODO: how to combine quadrature rules with the contours?
struct quadrature{shape}
    reference::shape
    nodes::Vector{Vector{Float64}}
    weights::Vector{Float64}
end

function T(z, D)
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

#  N + 1 : the number of quadrature nodes
#  ellip : the contour
#  D     : the size of the NEP
#  L     : a number between D and k                     
#  T     : the matrix of the NEP
function cont_int(ellip::ellipse, T::Function, D, L, N)

    A0 = zeros(ComplexF64, D, L)
    A1 = zeros(ComplexF64, D, L)

    # generate random matrix Vhat
    Vhat = randn(ComplexF64, D, L)

    # compute A0 and A1 (contour integrals with some quadrature rules)
	# here we use trapezoid rule on a ellipse
    delta = 2*pi/N
    for n = 0:N - 1
        theta = n * delta
        c = cos(theta)
        s = sin(theta)
		z1 = complex(ellip.semi_x * c, ellip.semi_y * s)
		dz = complex(ellip.semi_y * c, ellip.semi_x * s) / N 
		A0 = A0 + inv(T(ellip.center+z1, D)) * Vhat * dz
		A1 = A1 + z1 * inv(T(ellip.center+z1, D)) * Vhat * dz
    end

    # SVD of A0
    (V, Sig, W) = svd(A0)

    tol = 1.0e-8
    k = length(findall(abs.(Sig) .> tol))
    println("Find $(k) nonzero singular values.")

    # perform a rank test

    # compute the matrix B
	V_k = V[:, 1:k]
    Sig_k = Sig[1:k]
    W_k = W[:, 1:k]
    B = V_k' * A1 * W_k * diagm(1 ./ Sig_k)

    # solve the eigenvalue problem for B
    lambda = eigvals(B)
    V = eigvecs(B)
    lambda = lambda + ellip.center*ones(size(lambda))

    return lambda, V
end

ell1 = ellipse(complex(150.0, 0.0), 148.0, 148.0)

(Lambda, v) = cont_int(ell1, T, 400, 10, 50)
@show Lambda
