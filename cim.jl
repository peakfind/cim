#=---------------------------------------------------------
- cim.jl
-
-
---------------------------------------------------------=#

using LinearAlgebra

# contours that we plan to support
abstract type contour end

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


function cont_int(ellip::ellipse, D, L, N)

    A0 = zeros(ComplexF64, D, L)
    A1 = zeros(ComplexF64, D, L)

    # generate random matrix Vhat
    Vhat = randn(ComplexF64, D, L)

    # compute A0 and A1 (contour integrals with some quadrature rules)
    DeltaTheta = 2.0*pi/N
    for n = 0:N - 1
        theta = n*DeltaTheta
        c = cos(theta)
        s = sin(theta)
        
        z1 = ellip.semi_x*c + ellip.semi_y*s
    end




    # SVD of A0

    # perform a rank test

    # compute the matrix B

    # solve the eigenvalue problem for B
    
end