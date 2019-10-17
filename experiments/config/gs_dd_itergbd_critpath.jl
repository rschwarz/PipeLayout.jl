using JuMP
using SCIP
using GLPK

solver = GndStr.IterGBD(
    JuMP.with_optimizer(SCIP.Optimizer, display_width=139, limits_memory=5000.0),
    JuMP.with_optimizer(GLPK.Optimizer),
    addnogoods=false,
    addcritpath=true,
    maxiter=Int(1e9),
    timelimit=3600.0)
