using SCIP

solver = GndStr.IterTopo(
    SCIP.Optimizer(display_width=139, limits_memory=5000.0),
    SCIP.Optimizer(display_width=139, limits_memory=5000.0),
    maxiter=Int(1e9),
    timelimit=3600.0,
    debug=true)
