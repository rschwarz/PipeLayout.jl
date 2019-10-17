using JuMP
using SCIP

solver = GndStr.MINLP(
    JuMP.with_optimizer(SCIP.Optimizer, display_width=139, limits_memory=5000.0),
    timelimit=3600.0,
    writemodels=true)
