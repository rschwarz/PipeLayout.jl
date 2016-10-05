const account = "gas"

config    = ARGS[1]
instlist  = ARGS[2]
partition = ARGS[3]
results   = ARGS[4]

# create dir for results
outdir = joinpath(results, basename(instlist), basename(config))
mkpath(outdir)

"submit a job to SLURM"
function submit(key)
    account = "--account $account --partition $partition"
    limits = "--cpus-per-task 1"
    out, err = joinpath(outdir, "$key.log"), joinpath(outdir, "$key.err")
    output = "--output $out --error $err"
    options = "$account $limits $output"
    job = "julia run.jl $config $key"
    run(`sbatch $options $job`)
end

# one job per instance
open(instlist) do f
    for line in eachline(f)
        key = strip(line)
        submit(key)
    end
end
