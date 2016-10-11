immutable MINLP <: GroundStructureSolver
    debug::Bool
end

function make_minlp(inst::Instance, topo::Topology, solver::MINLP)
    m = Model()

    # TODO
end

function optimize(inst::Instance, topo::Topology, solver::MINLP)
    m = make_minlp(inst, topo, solver)
    status = solve(m)

    # TODO
end
