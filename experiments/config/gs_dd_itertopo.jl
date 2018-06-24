using SCIP

solver = GndStr.IterTopo(
    SCIPSolver("display/width", 139, "limits/memory", 5000.0),
    SCIPSolver("display/width", 139, "limits/memory", 5000.0),
    maxiter=Int(1e9),
    timelimit=3600.0,
    debug=true)
