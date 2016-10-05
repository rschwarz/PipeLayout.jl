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

println("solving $(instname) with $(config)")

tic()
println("vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv")
result = solve_with(instname, solver)
println("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
toc()

println("Status: $(result.status)")
