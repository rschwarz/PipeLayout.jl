using SCIP
using GLPK

solver = GndStr.IterGBD(
    SCIP.Optimizer(display_width=139, limits_memory=5000.0),
    GLPK.Optimizer(),
    addnogoods=true,
    addcritpath=false,
    maxiter=Int(1e9),
    timelimit=3600.0,
    debug=true)
