using SCIP
using GLPK

solver = GndStr.CallbackGBD(
    SCIPSolver(display_width=139, limits_memory=5000.0),
    GLPK.Optimizer(),
    addnogoods=true,
    addcritpath=true,
    timelimit=3600.0,
    debug=true)
