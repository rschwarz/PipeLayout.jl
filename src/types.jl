# Types for problem classes and solvers

export PipeLayoutSolver, GroundStructureSolver, JunctionLocationSolver
export PipeDimensioningSolver, optimize

"Base type for all solvers"
abstract type PipeLayoutSolver end

"Base type for all solvers depending on ground structure"
abstract type GroundStructureSolver <: PipeLayoutSolver end

"Base type for all solvers of junction location with fixed topology"
abstract type JunctionLocationSolver <: PipeLayoutSolver end

"Base type for all solvers of pipe dimensioning with fixed geometry"
abstract type PipeDimensioningSolver <: PipeLayoutSolver end

function optimize(inst::Instance, topo::Topology, solver)
    throw(ArgumentError("No method defined for solver type $(typeof(solver))"))
end
