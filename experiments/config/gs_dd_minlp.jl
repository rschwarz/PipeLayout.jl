using SCIP
solver = GndStr.MINLP(
    SCIPSolver("display/width", 139, "limits/memory", 5000),
    debug=true,
    timelimit=3600.0,
    writemodels=true)
