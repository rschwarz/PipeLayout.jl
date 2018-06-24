#! /usr/bin/env julia
using PipeLayout
using PipeLayout.GndStr

const ACCOUNT = "gas"

if length(ARGS) != 4
    println("usage: run_all.jl <config> <instlist> <partition> <results>")
    exit(1)
end

CONFIG    = abspath(ARGS[1])
INSTLIST  = ARGS[2]
PARTITION = ARGS[3]
RESULTS   = abspath(ARGS[4])

stem(path) = split(basename(path), ".")[1]

# location of executable
EXEPATH = abspath("run.jl")

# location of instances
INDIR = abspath(dirname(INSTLIST))

# create dir for results
OUTDIR = joinpath(RESULTS, stem(INSTLIST), stem(CONFIG))
mkpath(OUTDIR)

# conservative timelimit (in minutes)
include(CONFIG) # creates solver
@assert solver.timelimit < Inf
timelimit = round(Int, 2 * solver.timelimit / 60)

"submit a job to SLURM"
function submit(key)
    OUT = abspath("$key.log")
    options = ["--account=$ACCOUNT",
               "--partition=$PARTITION",
               "--cpus-per-task=1",
               "--time=$timelimit",
               "--signal=B:INT",
               "--output=$OUT",
               "--error=$OUT"]
    job = [EXEPATH, CONFIG, joinpath(INDIR, key)]
    run(`sbatch $options $job`)
end

# one job per instance, submit from result dir
open(INSTLIST) do f
    for line in eachline(f)
        key = strip(line)
        path = joinpath(OUTDIR, key)

        mkpath(path)
        cd(path)
        submit(key)
    end
end
