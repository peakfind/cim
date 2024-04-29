using Plots
# use plotlyjs() to generate interactive figure

# contours that we plan to support:
#     1. ellipse
#     2. circle
#     3. rectangle (not finished)
abstract type AbstractContour end

struct ellipse <: AbstractContour
    center::Vector{Float64} # TODO: store the vector{F64} or complex?
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

struct quadpts
    N::Int64                      # the number of the quadrature nodes
    nodes::AbstractArray{Float64} # quadrature nodes
end

function show_contr(contour::AbstractContour)
    if typeof(contour) == ellipse
        t = range(0, 2*pi, length=100)
        x1 = contour.center[1] .+ contour.semi_x*cos.(t)
        x2 = contour.center[2] .+ contour.semi_y*sin.(t)
        display(plot(x1, x2, linestyle=:dash, linecolor=:red, legend=:false, aspect_ratio=:equal))
    elseif typeof(contour) == circle
        t = range(0, 2*pi, length=100)
        x1 = contour.center[1] .+ contour.radius*cos.(t)
        x2 = contour.center[2] .+ contour.radius*sin.(t)
        display(plot(x1, x2, linestyle=:dash, linecolor=:red, legend=:false, aspect_ratio=:equal))
    # elseif typeof(contour) == rectangle
        # display(plot())
    else
        # throw some error : unsupported type of contour
        throw(ArgumentError("unspported type of contour"))
    end
end

# TODO: note that rectangle is different with ellipse and circle
# TODO: we may add the derivative of the parametrization in quadpts
function get_quadpts(contour::AbstractContour, num_quadpts::Int64)
    nodes = zeros(num_quadpts, 2)
    delta = 2*pi / num_quadpts

    if typeof(contour) == ellipse
        for i = 0:num_quadpts-1
            nodes[i + 1, 1] = contour.center[1] + contour.semi_x * cos(delta * i)
            nodes[i + 1, 2] = contour.center[2] + contour.semi_y * sin(delta * i)
        end
    elseif typeof(contour) == circle
        for i = 0:N-1
            nodes[i + 1, 1] = contour.center[1] + contour.radius * cos(delta * i) 
            nodes[i + 1, 2] = contour.center[1] + contour.radius * sin(delta * i)
        end
    else
        throw(ArgumentError("unspported type of contour"))
    end

    return quadpts(num_quadpts, nodes)
end

function show_quadpts(pts::quadpts)
    display(scatter!(pts.nodes[:,1], pts.nodes[:,2], markershape=:circle, markercolor=:blue))
end