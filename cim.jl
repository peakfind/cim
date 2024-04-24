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

#  N the number of quadrature nodes
function cont_int(ellip::ellipse, D, L, N)

    A0 = zeros(ComplexF64, D, L)
    A1 = zeros(ComplexF64, D, L)

    # generate random matrix Vhat
    V̂ = randn(ComplexF64, D, L)

    # compute A0 and A1 (contour integrals with some quadrature rules)
    DeltaTheta = 2.0*pi/N
    for n = 0:N - 1
        theta = n*DeltaTheta
        c = cos(theta)
        s = sin(theta)
        
        z1 = ellip.semi_x*c + ellip.semi_y*s
    end

    # SVD of A0
    (V, Σ, W) = svd(A0)

    tol = 1.0e-8
    k = length(findall(abs(Σ) .> tol))
    println("Find $(k) nonzero singular values.")

    # perform a rank test

    # compute the matrix B
    Vk = V[:, 1:k]
    Simgak = Σ[1:k]
    Wk = W[:, 1:k]
    B = Vk' * A1 * Wk * diagm(1 ./ Sigmak)

    # solve the eigenvalue problem for B
    λ = eigvals(B)
    V = eigvecs(B)

    return λ, V
    
end