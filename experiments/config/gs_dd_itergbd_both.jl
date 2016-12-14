using Gurobi: GurobiSolver

solver = GndStr.IterGBD(
    GurobiSolver(OutputFlag=0),
    GurobiSolver(OutputFlag=0),
    addnogoods=true,
    addcritpath=true,
    maxiter=Int(1e9),
    timelimit=3600.0,
    debug=true)
