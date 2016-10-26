# Types for problem classes and solvers

export PipeLayoutSolver, GroundStructureSolver, JunctionLocationSolver

"Base type for all solvers"
abstract PipeLayoutSolver

"Base type for all solvers depending on ground structure"
abstract GroundStructureSolver <: PipeLayoutSolver

"Base type for all solvers of junction location with fixed topology"
abstract JunctionLocationSolver <: PipeLayoutSolver
