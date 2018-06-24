using SCIP

solver = GndStr.CallbackTopo(
    SCIPSolver("display/width", 139, "limits/memory", 5000.0),
    SCIPSolver("display/width", 139, "limits/memory", 5000.0),
    timelimit=3600.0,
    debug=true)
