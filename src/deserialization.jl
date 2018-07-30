import JSON

# from https://github.com/JuliaLang/JSON.jl/issues/39
# but without `flatten` or `Nullable`
function deserialize(::Type{T}, json::AbstractString) where T
    deserialize(T, JSON.parse(json))
end

function deserialize(::Type{T}, json::Dict) where T
    T([deserialize(T, f, json) for f in fieldnames(T)]...)
end

function deserialize(::Type{Tf}, field::AbstractString, json::Dict) where Tf
    deserialize(Tf, json[field])
end

function deserialize(::Type{T}, field::Symbol, json::Dict) where T
    deserialize(fieldtype(T,field), string(field), json)
end

function deserialize(::Type{T}, i::Integer) where T<:Integer
    T(i)
end

function deserialize(::Type{T}, i::AbstractString) where T<:Integer
    parse(T, i)
end

function deserialize(::Type{T}, i::AbstractFloat) where T<:AbstractFloat
    T(i)
end

function deserialize(::Type{T}, i::AbstractString) where T<:AbstractFloat
    parse(T, i)
end

# my own extension to vectors
function deserialize(::Type{Vector{T}}, A::Vector{Any}) where T
    T[deserialize(T, a) for a in A]
end

function read_instance(prefix, key::String)
    open(joinpath(prefix, "$(key).instance.json")) do f
        PipeLayout.deserialize(Instance, readstring(f))
    end
end
read_instance(key::String) = read_instance(".", key)

function read_topology(prefix, key::String)
    open(joinpath(prefix, "$(key).topology.json")) do f
        PipeLayout.deserialize(Topology, readstring(f))
    end
end
read_topology(key::String) = read_topology(".", key)
