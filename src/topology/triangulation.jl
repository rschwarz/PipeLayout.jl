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

function triangulate(points::Matrix{Float64}, switches=TriangleSwitches())
    s = "c" # create mesh from point Cloud
    s *= convert(String, switches)
    tio = TriangulateIO(pointlist=points)
    result, _ =  Triangulate.triangulate(s, tio)
    return result
end

function triangulate(poly::TriangulateIO, switches=TriangleSwitches())
    s = "p" # create mesh from Polygon
    s *= convert(String, switches)
    result, _ = Triangulate.triangulate(s, poly)
    return result
end

function delaunay_triangulation(points::Matrix{Float64})
    s = TriangleSwitches(conforming_delaunay=true, maximum_steiner_points=0)
    return triangulate(points, s)
end

function refine(trimesh::TriangulateIO, switches=TriangleSwitches())
    # use edges as segments
    poly = TriangulateIO(pointlist=trimesh.pointlist,
                         segmentlist=trimesh.edgelist)
    return triangulate(poly, switches)
end

pointset_mean(array) = dropdims(mean(array, dims=2), dims=2)

abstract type TriangleCenter end
struct TriangleCentroid <: TriangleCenter end
struct TriangleIncenter <: TriangleCenter end # angle bisectors intersection
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

function edge_centers(trimesh::TriangulateIO)::Matrix{Float64}
    @unpack pointlist, edgelist = trimesh
    return pointset_mean(pointlist[:, edgelist])
end

"Compute all angles (in degrees) in a triangle."
function angles(pos::Matrix{Float64})::Vector{Float64}
    @assert size(pos) == (2, 3)
    a = norm(pos[:, 2] - pos[:, 3])
    b = norm(pos[:, 1] - pos[:, 3])
    c = norm(pos[:, 1] - pos[:, 2])
    α = acos((b^2 + c^2 - a^2) / (2*b*c))
    β = acos((a^2 + c^2 - b^2) / (2*a*c))
    γ = π - α - β
    return (360 / 2π) .* [α, β, γ]
end

"Compute triangle area, given points."
function area(pos::Matrix{Float64})::Float64
    @assert size(pos) == (2, 3)

    # compute side lengths from point coordinates
    a = norm(pos[:, 2] - pos[:, 3])
    b = norm(pos[:, 1] - pos[:, 3])
    c = norm(pos[:, 1] - pos[:, 2])

    # apply Heron's formula
    s = (a + b + c)/2
    return sqrt(s*(s - a)*(s - b)*(s - c))
end

"""Subdivide triangles in six parts, adding center and edge midpoints.

Skip triangles with angles smaller than given threshold.
"""
function refine_sixths(trimesh::TriangulateIO; minimum_angle=0.0, minimum_area=0.0,
                       center=TriangleCentroid())::TriangulateIO
    @unpack pointlist, edgelist, trianglelist = trimesh
    @assert size(pointlist, 1) == 2
    @assert size(trianglelist, 1) == 3
    @assert size(edgelist, 1) == 2

    # points are indexed in sequence
    #  - old points:       1:NP                  (num. points)
    #  - triangle centers: NP+1:NP+NST           (num. subdivided triangles)
    #  - edge midpoints:   NP+NST+1:NP+NST+NSE   (num. edges)
    NP = size(pointlist, 2)
    NT = size(trianglelist, 2)
    NE = size(edgelist, 2)

    # count subdivided triangles adjacent to edge
    counter = Dict{Tuple{Int, Int}, Int}()

    # triangle centers
    all_centers::Matrix{Float64} = triangle_centers(trimesh, center)
    center_points = Vector{Float64}[]
    center_edges = Vector{Int32}[]
    NST = 0
    for t in 1:NT
        triangle = trianglelist[:, t]
        positions = pointlist[:, triangle]
        if minimum(angles(positions)) < minimum_angle
            continue
        elseif area(positions) < minimum_area
            continue
        end

        NST += 1
        center_idx = NP + NST
        push!(center_points, all_centers[:, t])
        for corner in triangle
            push!(center_edges, [corner, center_idx])
        end

        # mark side edges for subdivision
        for (i, j) in ((1, 2), (2, 3), (3, 1))
            edge = minmax(triangle[i], triangle[j])
            counter[edge] = get!(counter, edge, 0) + 1
        end
    end

    # edge midpoints
    all_midpoints::Matrix{Float64} = edge_centers(trimesh)
    midpoints = Vector{Float64}[]
    midpoint_edges = Vector{Int32}[]

    NSE = 0
    for e in 1:NE
        # only subdivide edges that have subdivided triangles on both sides
        edge = minmax(edgelist[:, e]...)
        if get(counter, edge, 0) < 2
            # keep existing edge
            push!(midpoint_edges, edgelist[:, e])
        else
            # subdivide
            NSE += 1
            midpoint_idx = NP + NST + NSE
            push!(midpoints, all_midpoints[:, e])
            for point in edgelist[:, e]
                push!(midpoint_edges, [point, midpoint_idx])
            end
        end
    end

    # We are still missing the edges from the midpoints to the centers. But the
    # Triangulate library can take care of that...

    points = hcat(pointlist, center_points..., midpoints...)
    @assert size(points) == (2, NP + NST + NSE)
    edges = hcat(center_edges..., midpoint_edges...)
    @assert size(edges) == (2, 3NST + NE + NSE)
    poly = TriangulateIO(pointlist=points, segmentlist=edges)
    return refine(triangulate(poly))
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
