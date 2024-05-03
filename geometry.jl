using Plots
# use plotlyjs() to generate interactive figure

"""
contours that we plan to support:
     1. ellipse
     2. circle
     3. rectangle (not supported yet)
"""
abstract type AbstractContour end

"""
"""
struct ellipse <: AbstractContour
    center::Vector{Float64} 
    semi_x::Float64
    semi_y::Float64
end

"""
"""
struct rectangle <: AbstractContour
    center::Vector{Float64}
    semi_x::Float64
    semi_y::Float64
end

"""
"""
struct circle <: AbstractContour
    center::Vector{Float64}
    radius::Float64
end

"""
"""
struct quadpts
    N::Int64                     # the number of the quadrature nodes
    nodes::Matrix{Float64}       # quadrature nodes size of N x 2 
    nodes_prime::Matrix{Float64} # derivative of the parametrization size of N x 2
end

function show_contr(contour::ellipse)
    t = range(0, 2*pi, length=100)
    x1 = contour.center[1] .+ contour.semi_x*cos.(t)
    x2 = contour.center[2] .+ contour.semi_y*sin.(t)
    display(plot(x1, x2, linestyle=:dash, linecolor=:red, legend=:false, aspect_ratio=:equal))
end

function show_contr(contour::circle)
    t = range(0, 2*pi, length=100)
    x1 = contour.center[1] .+ contour.radius*cos.(t)
    x2 = contour.center[2] .+ contour.radius*sin.(t)
    display(plot(x1, x2, linestyle=:dash, linecolor=:red, legend=:false, aspect_ratio=:equal))
end

"""
    get_quadpts(contour::ellipse or circle, num_quadpts::Int64)
    -------------------------------------------------------
    Input:
        contour    : the contour that we discretize 
        num_quadpts: the number of the quadrature nodes
    Output:
        quadpts    : the struct quadpts which contains the 
                     information of the quadrature nodes

    TODO: note that rectangle is different with ellipse and circle
"""
function get_quadpts(contour::ellipse, num_quadpts::Int64)
    nodes = zeros(num_quadpts, 2)
    nodes_prime = zeros(num_quadpts, 2)
    delta = 2*pi / num_quadpts

    for i = 0:num_quadpts-1
        nodes[i + 1, 1] = contour.center[1] + contour.semi_x * cos(delta * i)
        nodes[i + 1, 2] = contour.center[2] + contour.semi_y * sin(delta * i)
        nodes_prime[i + 1, 1] = -contour.semi_x * sin(delta * i)
        nodes_prime[i + 1, 2] = contour.semi_y * cos(delta * i)
    end

    return quadpts(num_quadpts, nodes, nodes_prime)
end

function get_quadpts(contour::circle, num_quadpts::Int64)
    nodes = zeros(num_quadpts, 2)
    nodes_prime = zeros(num_quadpts, 2)
    delta = 2*pi / num_quadpts

    for i = 0:num_quadpts-1
        nodes[i + 1, 1] = contour.center[1] + contour.radius * cos(delta * i) 
        nodes[i + 1, 2] = contour.center[1] + contour.radius * sin(delta * i)
        nodes_prime[i + 1, 1] = -contour.radius * sin(delta * i)
        nodes_prime[i + 1, 2] = contour.radius * cos(delta * i)
    end

    return quadpts(num_quadpts, nodes, nodes_prime)
end

function show_quadpts(pts::quadpts)
    display(scatter!(pts.nodes[:,1], pts.nodes[:,2], markershape=:circle, markercolor=:gray, markersize=2, aspect_ratio=:equal))
end

function show_eigs(eigs::AbstractArray)
    display(scatter!(real(eigs), imag(eigs), markershape=:xcross, markercolor=:black, markersize=3, aspect_ratio=:equal))
end