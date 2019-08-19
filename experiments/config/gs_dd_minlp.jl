using SCIP
solver = GndStr.MINLP(
    SCIP.Optimizer(display_width=139, limits_memory=5000.0),
    debug=true,
    timelimit=3600.0,
    writemodels=true)
