using SCIP
using GLPKMathProgInterface

solver = GndStr.CallbackGBD(
    SCIPSolver("display/width", 139, "limits/memory", 5000.0),
    GLPKSolverLP(),
    addnogoods=true,
    addcritpath=true,
    timelimit=3600.0,
    debug=true)
