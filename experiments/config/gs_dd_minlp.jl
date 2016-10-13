using SCIP
solver = MINLP(SCIPSolver(),
               debug=true,
               timelimit=3600.0)
