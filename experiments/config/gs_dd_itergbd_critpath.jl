using SCIP
using GLPKMathProgInterface

solver = GndStr.IterGBD(
    SCIPSolver("display/width", 139, "limits/memory", 5000.0),
    GLPKSolverLP(),
    addnogoods=false,
    addcritpath=true,
    maxiter=Int(1e9),
    timelimit=3600.0,
    debug=true)
