include("../cim.jl")
# we may need high order contour integral method to check this example

function F(z::Complex)
    return [exp(im * z * z) 1.0; 1.0 1.0]
end

n = 200 # nbr of quadrature nodes
circ = circle([0.0, 0.0], 3.0)
pts = get_quadpts(circ, n)
lambda = contr_int(pts, F, 2, 1)

@show lambda