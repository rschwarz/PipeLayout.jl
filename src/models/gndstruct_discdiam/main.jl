using MathOptInterface
using JuMP
using SparseArrays

const MOI = MathOptInterface

include("common.jl")

# iterative algorithms (decomposition with two separate models)
include("itergbd.jl")
include("itertopo.jl")

# all-in-one integrated model
include("minlp.jl")

# callback-based algorithms: keep one master problem
# include("cbtopo.jl")  # needs itertopo.jl
# include("cbgbd.jl")   # needs itergbd.jl
