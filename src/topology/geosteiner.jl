"Check whether geosteiner executables are in PATH"
function has_geosteiner()
    for cmd in [`efst -h`, `bb -h`]
        try
            run(pipeline(cmd, stdout=devnull, stderr=devnull))
        catch exc
            isa(exc, Base.IOError) && return false
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

    tnodes = copy(nodes)
    tarcs = Arc[]

    output, input, process = readandwrite(pipeline(`efst`, `bb`))
    for node in nodes
        write(input, "$(node.x) $(node.y)\n")
    end
    close(input)

    # extract relevant lines from postscript output of GeoSteiner
    line_tokens = Array{String}[]
    _started = false
    for line in eachline(output)
        # are we in the good part yet?
        if !_started && !contains(line, " % fs")
            continue
        else
            _started = true
        end

        # are we still in the good part?
        if contains(line, "Euclidean SMT")
            break
        end

        # skip the noise
        if contains(line, " % fs")
            continue
        end

        # extract the values
        push!(line_tokens, split(line))
    end

    for tokens in line_tokens
        @assert length(tokens) == 5
        @assert tokens[5] == "S"
        tail, head = 0, 0
        if tokens[2] == "T"
            tail = parse(Int, tokens[1]) + 1
        else
            stein = Node(parse(Float64, tokens[1]), parse(Float64, tokens[2]))
            tail = findfirst(tnodes, stein)
            if tail == 0
                push!(tnodes, stein)
                tail = length(tnodes)
            end
        end
        if tokens[4] == "T"
            head = parse(Int, tokens[3]) + 1
        else
            stein = Node(parse(Float64, tokens[3]), parse(Float64, tokens[4]))
            head = findfirst(tnodes, stein)
            if head == 0
                push!(tnodes, stein)
                head = length(tnodes)
            end
        end

        push!(tarcs, Arc(tail,head))
    end

    @assert(length(tnodes) <= 2*length(nodes) - 2)
    @assert(length(tarcs) == length(tnodes) - 1)

    Topology(tnodes, tarcs)
end
