using JuMP
using SCIP
using GLPK

solver = GndStr.CallbackGBD(
    JuMP.with_optimizer(SCIP.Optimizer, display_width=139, limits_memory=5000.0),
    JuMP.with_optimizer(GLPK.Optimizer),
    addnogoods=true,
    addcritpath=true,
    timelimit=3600.0)
