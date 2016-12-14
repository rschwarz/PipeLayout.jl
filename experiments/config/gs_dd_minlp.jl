using SCIP
solver = GndStr.MINLP(
    SCIPSolver(),
    debug=true,
    timelimit=3600.0)
