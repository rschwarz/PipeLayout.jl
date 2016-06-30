using JuMP

# TODO: actually build the model, think about naming
function build_pd_model(inst::Instance, topo::Topology)
    n = length(instance.nodes)
    length(topo.nodes) == n || throw(ArgumentError("Steiner nodes not allowed"))
    d = length(instance.diameters)

    q = flow(inst, topo)

    model = Model()

    # squared pressure variables at nodes
    lb, ub = [b.lb for b in inst.bounds], [b.ub for b in inst.bounds]
    @variable(model, lb <= π[1:n] <= ub)

    # relative length of pipe segments with specific diameter
    @variable(model, 0 <= l[])


end
