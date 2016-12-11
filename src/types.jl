# Types for problem classes and solvers

export PipeLayoutSolver, GroundStructureSolver, JunctionLocationSolver
export PipeDimensioningSolver, optimize

"Base type for all solvers"
abstract PipeLayoutSolver

"Base type for all solvers depending on ground structure"
abstract GroundStructureSolver <: PipeLayoutSolver

"Base type for all solvers of junction location with fixed topology"
abstract JunctionLocationSolver <: PipeLayoutSolver

"Base type for all solvers of pipe dimensioning with fixed geometry"
abstract PipeDimensioningSolver <: PipeLayoutSolver

function optimize(inst::Instance, topo::Topology, solver)
    throw(ArgumentError("No method defined for solver type $(typeof(solver))"))
end
