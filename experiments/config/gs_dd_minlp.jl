using JuMP
using SCIP

solver = GndStr.MINLP(
    JuMP.optimizer_with_attributes(
        SCIP.Optimizer,
        "display/width" => 139,
        "limits/memory" => 5000.0
    ),
    timelimit=3600.0,
    writemodels=true)
