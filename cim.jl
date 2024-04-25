#=---------------------------------------------------------
- cim.jl
- TODO: we should write a script to implement the method first
-
---------------------------------------------------------=#

using LinearAlgebra

# contours that we plan to support
abstract type contour end

# TODO: Can we provide a generic method for different types of contours?
struct ellipse <: contour
    center::Vector{Float64}
    semi_x::Float64
    semi_y::Float64
end

struct rectangle <: contour
    center::Vector{Float64}
    semi_x::Float64
    semi_y::Float64
end

struct circle <: contour
    center::Vector{Float64}
    radius::Float64
end

# TODO: how to combine quadrature rules with the contours?
struct quadrature{shape}
    reference::shape
    nodes::Vector{Vector{Float64}}
    weights::Vector{Float64}
end

#  N     : the number of quadrature nodes
#  ellip : the contour
#  D     :
#  L     :
#  T     : the matrix of the NEP
function cont_int(ellip::ellipse, T, D, L, N)

    A0 = zeros(ComplexF64, D, L)
    A1 = zeros(ComplexF64, D, L)

    # generate random matrix Vhat
    Vhat = randn(ComplexF64, D, L)

    # compute A0 and A1 (contour integrals with some quadrature rules)
	# here we use trapezoid rule on a ellipse
    delta = 2.0*pi/N
    for n = 0:N - 1
        theta = n * delta
        c = cos(theta)
        s = sin(theta)
		z1 = complex(ellip.semi_x * c, ellip.semi_y * s)
		dz = complex(ellip.semi_y * c, ellip.semi_x * s) / N 
		A0 = A0 + T * Vhat * dz
		A1 = A1 + z1 * T * Vhat * dz
    end

    # SVD of A0
    (V, Sig, W) = svd(A0)

    tol = 1.0e-8
    k = length(findall(abs(Sig) .> tol))
    println("Find $(k) nonzero singular values.")

    # perform a rank test

    # compute the matrix B
    Vk = V[:, 1:k]
    Sigk = Sig[1:k]
    Wk = W[:, 1:k]
    B = Vk' * A1 * Wk * diagm(1 ./ Sigk)

    # solve the eigenvalue problem for B
    lambda = eigvals(B)
    V = eigvecs(B)

    return lambda, V
end
