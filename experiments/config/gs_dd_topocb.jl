using JuMP
using SCIP

solver = GndStr.CallbackTopo(
    JuMP.with_optimizer(SCIP.Optimizer, display_width=139, limits_memory=5000.0),
    JuMP.with_optimizer(SCIP.Optimizer, display_width=139, limits_memory=5000.0),
    timelimit=3600.0,
    debug=true)
