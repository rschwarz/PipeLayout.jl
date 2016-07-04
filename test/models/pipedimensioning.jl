import PipeLayout: pipedim_model
using JuMP

facts("check dimensioning on single pipe") do
    inst, topo = single_pipe()
    model, π, l = pipedim_model(inst, topo)
    status = solve(model)
    @show status

    psol = getvalue(π).^0.5
    @fact psol[1] --> roughly_within(60.0, 80.0)
    @fact psol[2] --> roughly_within(60.0, 80.0)
    @fact psol[1] --> greater_than(psol[2])

    lsol = getvalue(l)
    @fact size(lsol) --> (1,3) # one pipe with 3 segments
    @fact sum(lsol[1,:]) --> roughly(1.0)
    @fact lsol[1,1] --> roughly_within(0.0, 1.0)
    @fact lsol[1,2] --> roughly_within(0.0, 1.0)
    @fact lsol[1,3] --> roughly_within(0.0, 1.0)
    @fact lsol[1,:] --> is_SOS2
end
