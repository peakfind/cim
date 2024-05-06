#=-------------------------------------------------------------------------------------
- cim.jl
- 
- references:
- 1. Wolf-Jurgen Beyn, An integral method for solving NEPs, Linear Algebra Appl., 2012.
- 2. Stefan Guttel, Francoise Tisseur, The NEP, Acta Numerica, 2017.
-
-------------------------------------------------------------------------------------=#

using LinearAlgebra
include("geometry.jl")

"""
NEP: the matrix function of the nonlinear eigenvalue problem
D: the size of the NEP
l: a number between k(the number of eigenvalues inside the contour) and D, k <= L <= D
"""
function contr_int(pts::quadpts, NEP::Function, D, l::Int64)
    # Preallocate arrays
    A0 = zeros(ComplexF64, D, l)
    A1 = zeros(ComplexF64, D, l)
    Vhat = randn(ComplexF64, D, l)

    # Compute A0 and A1 with trapezoid rule
    for j in 1:pts.N
        z = complex(pts.nodes[j, 1], pts.nodes[j, 2])
        z_prime = complex(pts.nodes_prime[j, 1], pts.nodes_prime[j, 2])
        invNEP_Vhat = NEP(z) \ Vhat
        A0 .+= invNEP_Vhat * z_prime
        A1 .+= A1 + invNEP_Vhat * z * z_prime
    end
    A0 ./= (pts.N*im)
    A1 ./= (pts.N*im)

    # Compute the SVD of A0
    (V, Sigma, W) = svd(A0)

    # Determine the number of nonzero singular values 
    tol = 1.0e-12
    k = count(Sigma ./ Sigma[1] .> tol)

    # Compute the matrix B 
    Vk = V[:,1:k]
    Sigk = Sigma[1:k]
    Wk = W[:,1:k]

    # Diagonal is more efficient
    B = (Vk' * A1 * Wk) * Diagonal(1 ./ Sigk)

    # Compute the eigenvalues of B 
    lambda = eigvals(B)

    return lambda
end

"""
contour integral method with high-order moments 
pbar: the number of the moments (for p = 0, ..., pbar)

In this function, we follow the notations in Stefan Guttel, Francoise Tisseur, Acta Numerica, 2017
r,l : size of the probing matrices
"""
function contr_int_ho(pts::quadpts, NEP::Function, D::Int64, l::Int64, r::Int64, pbar::Int64)
    # initialization
    L = randn(ComplexF64, D, l)
    R = randn(ComplexF64, D, r)
    A = zeros(ComplexF64, l, r, 2*pbar)# contains all the moments p = 0, ..., 2*pbar - 1

    # compute moments
    # TODO: we need use NEP\R not the inv()
    for j in 1:pts.N
        z = complex(pts.nodes[j, 1], pts.nodes[j, 2])
        z_prime = complex(pts.nodes_prime[j, 1], pts.nodes_prime[j, 2])
        for p in 1:2*pbar
            A[:,:,p] = A[:,:,p] + (L' * (NEP(z)\R)) * (z^p) * z_prime
        end
    end

    A = A / (pts.N * im)

    # compute B0 and B1
    B0 = zeros(ComplexF64, l*pbar, r*pbar)
    B1 = zeros(ComplexF64, l*pbar, r*pbar)

    for j = 1:pbar
        B0[1+(j-1)*l:j*l,:] = reshape(A[:,:,j:j+pbar-1], l,:)
        B1[1+(j-1)*l:j*l,:] = reshape(A[:,:,j+1:j+pbar], l,:) 
    end

    # svd of B0
    println("Compute the SVD of B0.")
    (V, Sigma, W) = svd(B0)

    # mbar: the number of the eigs inside the contour
    tol = 1.0e-12
    mbar = count(Sigma/Sigma[1] .> tol)
    println("Find $(mbar) nonzero singular values.")

    # compute the matrix M
    Vmbar = V[:,1:mbar]
    Sigmbar = Sigma[1:mbar]
    Wmbar = W[:,1:mbar]

    M = (Vmbar' * B1 * Wmbar) * Diagonal(1 ./ Sigmbar)

    # compute the eigenvalues of M
    lambda = eigvals(M)

    return lambda
end

function contr_int_ho_gpt(pts::quadpts, NEP::Function, D::Int64, l::Int64, r::Int64, pbar::Int64)
    # Preallocate arrays
    L = randn(ComplexF64, D, l)
    R = randn(ComplexF64, D, r)
    A = zeros(ComplexF64, l, r, 2*pbar)# contains all the moments p = 0, ..., 2*pbar - 1

    # Compute moments (We need use NEP\R not the inv())
    for j in 1:pts.N
        z = complex(pts.nodes[j, 1], pts.nodes[j, 2])
        z_prime = complex(pts.nodes_prime[j, 1], pts.nodes_prime[j, 2])
        L_invNEP_R = L' * (NEP(z) \ R)
        for p in 1:2*pbar
            A[:,:,p] .+= L_invNEP_R * (z^p) * z_prime
        end
    end

    A ./= (pts.N * im)

    # Compute B0 and B1
    B0 = zeros(ComplexF64, l*pbar, r*pbar)
    B1 = zeros(ComplexF64, l*pbar, r*pbar)

    for j = 1:pbar
        B0[1+(j-1)*l:j*l,:] = reshape(A[:,:,j:j+pbar-1], l,:)
        B1[1+(j-1)*l:j*l,:] = reshape(A[:,:,j+1:j+pbar], l,:) 
    end

    # Compute the SVD of B0
    (V, Sigma, W) = svd(B0)

    # Determine the number of nonzero singular values
    # mbar: the number of the eigs inside the contour
    tol = 1.0e-12
    mbar = count(Sigma ./ Sigma[1] .> tol)

    # Compute the matrix M
    Vmbar = V[:,1:mbar]
    Sigmbar = Sigma[1:mbar]
    Wmbar = W[:,1:mbar]
    M = (Vmbar' * B1 * Wmbar) * Diagonal(1 ./ Sigmbar)

    # Compute the eigenvalues of M
    lambda = eigvals(M)

    return lambda
end