using SCIP
solver = GndStr.MINLP(
    SCIPSolver("display/width", 139),
    debug=true,
    timelimit=3600.0,
    writemodels=true)
