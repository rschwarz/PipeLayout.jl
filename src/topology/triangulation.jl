using LinearAlgebra
using Statistics

using LightGraphs
using Parameters
using Triangulate

@with_kw struct TriangleSwitches
    minimum_angle = nothing           # numeric
    maximum_area = nothing            # numeric
    conforming_delaunay = false
    maximum_steiner_points = nothing  # numeric
    convex_hull = false
    no_boundary_steiner_points = false
    check_consistency = true
    quiet = true
end

function Base.convert(::Type{String}, x::TriangleSwitches)
    s = ""
    s *= "e" # always output edges
    if x.minimum_angle !== nothing
        s *= "q$(x.minimum_angle)"
    end
    if x.maximum_area !== nothing
        s *= "a$(x.maximum_area)"
    end
    if x.conforming_delaunay
        s *= "D"
    end
    if x.maximum_steiner_points !== nothing
        s *= "S$(x.maximum_steiner_points)"
    end
    if x.convex_hull
        s *= "c"
    end
    if x.no_boundary_steiner_points
        s *= "Y"
    end
    if x.check_consistency
        s *= "C"
    end
    if x.quiet
        s *= "Q"
    end
    return s
end

function Base.convert(::Type{Topology}, x::TriangulateIO)
    nodes = Node.(eachcol(x.pointlist))
    arcs = Arc.(eachcol(x.edgelist))
    return Topology(nodes, arcs)
end

function triangulate(points::Matrix{Float64}, switches::TriangleSwitches)
    s = "c" # create mesh from point Cloud
    s *= convert(String, switches)
    tio = TriangulateIO(pointlist=points)
    result, _ =  Triangulate.triangulate(s, tio)
    return result
end

function triangulate(poly::TriangulateIO, switches::TriangleSwitches)
    s = "p" # create mesh from Polygon
    s *= convert(String, switches)
    result, _ = Triangulate.triangulate(s, poly)
    return result
end

function delaunay_triangulation(points::Matrix{Float64})
    s = TriangleSwitches(conforming_delaunay=true, maximum_steiner_points=0)
    return triangulate(points, s)
end

pointset_mean(array) = dropdims(mean(array, dims=2), dims=2)

abstract type TriangleCenter end
struct TriangleCentroid <: TriangleCenter end
struct TriangleIncenter <: TriangleCenter end
struct TriangleCircumcenter <: TriangleCenter end # yield Voronoi points!

function triangle_centers(trimesh::TriangulateIO, ::TriangleCentroid)
    @unpack pointlist, trianglelist = trimesh
    return pointset_mean(pointlist[:, trianglelist])
end

function triangle_centers(trimesh::TriangulateIO, ::TriangleIncenter)
    @unpack pointlist, trianglelist = trimesh
    centers = []
    for t in eachcol(trianglelist)
        corners = pointlist[:, t]

        a = norm(corners[:, 2] - corners[:, 3])
        b = norm(corners[:, 1] - corners[:, 3])
        c = norm(corners[:, 1] - corners[:, 2])
        # based on barycentric coordinates a:b:c
        @show [a; b; c] corners
        incenter = corners * [a, b, c] ./ (a + b + c)
        push!(centers, incenter)
    end

    return hcat(centers...)
end

function triangle_centers(trimesh::TriangulateIO, ::TriangleCircumcenter)
    @unpack pointlist, trianglelist = trimesh
    centers = []
    for t in eachcol(trianglelist)
        corners = pointlist[:, t]
        Ax, Ay = corners[:, 1]
        Bx, By = corners[:, 2]
        Cx, Cy = corners[:, 3]
        D = 2 * ( Ax * (By - Cy) + Bx * (Cy - Ay) + Cx * (Ay - By) )
        Ux = ((Ax^2 + Ay^2)*(By - Cy)) + ((Bx^2 + By^2)*(Cy - Ay)) + ((Cx^2 + Cy^2)*(Ay - By))
        Uy = ((Ax^2 + Ay^2)*(Cx - Bx)) + ((Bx^2 + By^2)*(Ax - Cx)) + ((Cx^2 + Cy^2)*(Bx - Ax))
        circumcenter = [Ux, Uy] ./ D
        push!(centers, circumcenter)
    end

    return hcat(centers...)
end

function antiparallel_digraph(trimesh::TriangulateIO)
    @unpack pointlist, edgelist = trimesh
    graph = SimpleDiGraph(size(pointlist, 2))
    for (s, t) in eachcol(edgelist)
        add_edge!(graph, s, t)
        add_edge!(graph, t, s)
    end
    return graph
end
