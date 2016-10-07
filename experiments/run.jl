#! /usr/bin/env julia

# raise exception when receiving SIGINT
ccall(:jl_exit_on_sigint, Void, (Cint,), 0)

using PipeLayout
using PipeLayout.GndStructDiscDiam

function solve_with(instance::AbstractString, solver::GroundStructureSolver)
    inst = PipeLayout.read_instance(instance)
    topo = PipeLayout.read_topology(instance)

    optimize(inst, topo, solver)
end

config = ARGS[1]
instname = ARGS[2]

include(config) # creates solver

println("solving $instname with $(basename(config))")

tic()
println("vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv")
try
    result = solve_with(instname, solver)
catch ex
    if isa(ex, InterruptException)
        println("-- was interrupted --")
    else
        rethrow()
    end
end
println("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
toc()

println("Status: $(result.status)")
