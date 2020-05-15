using JuMP
using SCIP
using GLPK

solver = GndStr.CallbackGBD(
    JuMP.optimizer_with_attributes(
        SCIP.Optimizer,
        "display/width" => 139,
        "limits/memory" => 5000.0
    ),
    GLPK.Optimizer,
    addnogoods=true,
    addcritpath=true,
    timelimit=3600.0)
