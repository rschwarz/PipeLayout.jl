"Check whether geosteiner executables are in PATH"
function has_geosteiner()
    for cmd in [`efst -h`, `bb -h`]
        try
            run(pipeline(cmd, DevNull))
        catch exc
            isa(exc, Base.UVError) && return false
        end
    end
    true
end


"""
Compute a Euclidean Steiner Minimal Tree for given nodes.

Is implemented by "shelling out" to GeoSteiner.
"""
function euclidean_steiner_tree(nodes::Vector{Node})
    @assert has_geosteiner() "Need GeoSteiner executables (efts, bb) in PATH!"

    output, input, process = readandwrite(pipeline(`efst`, `bb`))
    for node in nodes
        write(input, "$(node.x) $(node.y)\n")
    end
    @show input, output, process
    close(input)
    @show input, output, process
    println(readall(output))
    @show input, output, process
end
