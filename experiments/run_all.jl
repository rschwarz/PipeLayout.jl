#! /usr/bin/env julia

const ACCOUNT = "gas"

CONFIG    = ARGS[1]
INSTLIST  = ARGS[2]
PARTITION = ARGS[3]
RESULTS   = ARGS[4]

stem(path) = split(basename(path), ".")[1]

# location of instances
INDIR = abspath(dirname(INSTLIST))

# create dir for results
OUTDIR = joinpath(RESULTS, stem(INSTLIST), stem(CONFIG))
mkpath(OUTDIR)

"submit a job to SLURM"
function submit(key)
    OUT = joinpath(OUTDIR, "$key.log")
    ERR = joinpath(OUTDIR, "$key.err")
    options = ["--account=$ACCOUNT",
               "--partition=$PARTITION",
               "--cpus-per-task=1",
               "--time=2:05:00",
               "--signal=SIGINT",
               "--output=$OUT",
               "--error=$ERR"]
    job = ["run.jl", abspath(CONFIG), joinpath(INDIR, key)]
    run(`sbatch $options $job`)
end

# one job per instance
open(INSTLIST) do f
    for line in eachline(f)
        key = strip(line)
        submit(key)
    end
end
