import JSON

# from https://github.com/JuliaLang/JSON.jl/issues/39
# but without `flatten` or `Nullable`
function deserialize{T}(::Type{T}, json::AbstractString)
    deserialize(T, JSON.parse(json))
end

function deserialize{T}(::Type{T}, json::Dict)
    T([deserialize(T, f, json) for f in fieldnames(T)]...)
end

function deserialize{Tf}(::Type{Tf}, field::AbstractString, json::Dict)
    deserialize(Tf, json[field])
end

function deserialize{T}(::Type{T}, field::Symbol, json::Dict)
    deserialize(fieldtype(T,field), string(field), json)
end

function deserialize{T<:Integer}(::Type{T}, i::Integer)
    T(i)
end

function deserialize{T<:Integer}(::Type{T}, i::AbstractString)
    parse(T, i)
end

function deserialize{T<:AbstractFloat}(::Type{T}, i::AbstractFloat)
    T(i)
end

function deserialize{T<:AbstractFloat}(::Type{T}, i::AbstractString)
    parse(T, i)
end

# my own extension to vectors
function deserialize{T}(::Type{Vector{T}}, A::Vector{Any})
    T[deserialize(T, a) for a in A]
end

"convenience to load json files with instance and topology"
function read_files(prefix, key::String)
    res = []
    open(joinpath(prefix, "$(key).instance.json")) do f
        push!(res, PipeLayout.deserialize(Instance, readstring(f)))
    end
    open(joinpath(prefix, "$(key).topology.json")) do f
        push!(res, PipeLayout.deserialize(Topology, readstring(f)))
    end
    res
end
