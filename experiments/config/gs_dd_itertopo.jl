using Gurobi: GurobiSolver

solver = IterTopo(GurobiSolver(OutputFlag=0),
                  GurobiSolver(OutputFlag=0),
                  maxiter=Int(1e9),
                  timelimit=3600.0,
                  debug=true)
