using JuMP

# TODO: actually build the model, think about naming
function build_pd_model(inst::Instance, topo::Topology)
    q = flow(inst, topo)
    m = Model()

end
