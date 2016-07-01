import PipeLayout: pipedim_model
using JuMP

facts("check dimensioning on single pipe") do
    inst, topo = single_pipe()
    model, π, l = pipedim_model(inst, topo)
    status = solve(model)
    @show status

    # TODO: do proper tests
end
