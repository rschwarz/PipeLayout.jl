# Types for problem classes and solvers

export PipeLayoutSolver, GroundStructureSolver

"Base type for all solvers"
abstract PipeLayoutSolver

"Base type for all solvers depending on ground structure"
abstract GroundStructureSolver <: PipeLayoutSolver
