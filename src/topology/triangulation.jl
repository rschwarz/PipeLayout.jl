using LightGraphs
using Parameters
using TriangleMesh

@with_kw struct TriangleSwitches
    minimum_angle = nothing           # numeric
    maximum_area = nothing            # numeric
    conforming_delaunay = false
    maximum_steiner_points = nothing  # numeric
    no_boundary_steiner_points = false
    check_consistency = true
    quiet = true
end

function Base.convert(String, x::TriangleSwitches)
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

"Create polygon from points and edges (used as segments)."
function make_pslg(points, edges)
    poly = Polygon_pslg(size(points, 1), 0, 0, size(edges, 1), 0)
    set_polygon_point!(poly, points)
    set_polygon_segment!(poly, edges)
    return poly
end

make_pslg(trimesh::TriMesh) = make_pslg(trimesh.point', trimesh.edge')

function triangulate(points::Matrix{Float64}, switches::TriangleSwitches)::TriMesh
    s = "c" # create mesh from point Cloud
    s *= convert(String, switches)
    return create_mesh(points, s)
end

function triangulate(poly::Polygon_pslg, switches::TriangleSwitches)::TriMesh
    s = "p" # create mesh from Polygon
    s *= convert(String, switches)
    return create_mesh(poly, s)
end

function delaunay_triangulation(points::Matrix{Float64})::TriMesh
    switches = TriangleSwitches(conforming_delaunay=true, maximum_steiner_points=0)
    return triangulate(points, switches)
end

pointset_mean(array) = dropdims(mean(array, dims=2), dims=2)

abstract type TriangleCenter end
struct TriangleCentroid <: TriangleCenter end
struct TriangleIncenter <: TriangleCenter end
struct TriangleCircumcenter <: TriangleCenter end # yield Voronoi points!

function triangle_centers(trimesh, ::TriangleCentroid)
    @unpack point, cell = trimesh
    return pointset_mean(point'[cell', :])
end

function triangle_centers(trimesh, ::TriangleIncenter)
    @unpack point, cell = trimesh
    centers = []
    for t in eachrow(cell')
        corners = point'[t, :]

        a = norm(corners[2, :] - corners[3, :])
        b = norm(corners[1, :] - corners[3, :])
        c = norm(corners[1, :] - corners[2, :])
        # based on barycentric coordinates a:b:c
        incenter = [a b c] * corners ./ (a + b + c)
        push!(centers, incenter)
    end

    return vcat(centers...)
end

function triangle_centers(trimesh, ::TriangleCircumcenter)
    @unpack point, cell = trimesh
    centers = []
    for t in eachrow(cell')
        corners = point'[t, :]
        Ax, Ay = corners[1, :]
        Bx, By = corners[2, :]
        Cx, Cy = corners[3, :]
        D = 2 * ( Ax * (By - Cy) + Bx * (Cy - Ay) + Cx * (Ay - By) )
        Ux = ((Ax^2 + Ay^2)*(By - Cy)) + ((Bx^2 + By^2)*(Cy - Ay)) + ((Cx^2 + Cy^2)*(Ay - By))
        Uy = ((Ax^2 + Ay^2)*(Cx - Bx)) + ((Bx^2 + By^2)*(Ax - Cx)) + ((Cx^2 + Cy^2)*(Bx - Ax))
        circumcenter = [Ux Uy] ./ D
        push!(centers, circumcenter)
    end

    return vcat(centers...)
end

function antiparallel_digraph(trimesh)
    points, edges = trimesh.point', trimesh.edge'
    graph = SimpleDiGraph(size(points, 1))
    for e in 1:size(edges, 1)
        s, t = edges[e, :]
        add_edge!(graph, s, t)
        add_edge!(graph, t, s)
    end
    return graph
end
