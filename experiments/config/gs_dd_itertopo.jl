using JuMP
using SCIP

solver = GndStr.IterTopo(
    JuMP.optimizer_with_attributes(
        SCIP.Optimizer,
        "display/width" => 139,
        "limits/memory" => 5000.0
    ),
    JuMP.optimizer_with_attributes(
        SCIP.Optimizer,
        "display/width" => 139,
        "limits/memory" => 5000.0
    ),
    maxiter=Int(1e9),
    timelimit=3600.0)
