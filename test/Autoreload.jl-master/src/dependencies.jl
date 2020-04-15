function topological_sort{T}(outgoing::Dict{T, Vector{T}})
    # kahn topological sort
    # todo: use better data structures
    order = T[]
    incoming = Dict{T,Vector}()
    for (n, children) in outgoing
        for child in children
            if haskey(incoming, child)
                push!(incoming[child], n)
            else
                incoming[child] = [n]
            end
        end
        if !haskey(incoming, n)
            incoming[n] = T[]
        end
    end
    S = T[]

    outgoing_counts = Dict{T,Int}()
    for (child, parents) in outgoing
        if isempty(parents)
            push!(S, child)
        end
        outgoing_counts[child] = length(parents)
    end
    while !isempty(S)
        n = pop!(S)
        push!(order, n)
        children = incoming[n]
        for child in children
            outgoing_counts[child] -= 1
            if outgoing_counts[child] == 0
                push!(S, child)
            end
        end
    end
    if length(order) != length(outgoing)
        error("Cyclic dependencies detected")
    end
    return order
end


function get_dependency_graph()
    deps = [filename=>afile.deps for (filename, afile) in files]
    return deps
end

function _extract_deps(e::Expr, deps::Vector{Dep})
    if e.head == :call
        if e.args[1] == :include
            if isa(e.args[2], AbstractString)
                push!(deps, Dep(false, e.args[2]))
            end
        # elseif e.args[1] == :require
        #     if isa(e.args[2], AbstractString)
        #         push!(deps, Dep(true, e.args[2]))
        #     end
        end
    elseif e.head == :import
        # push!(deps, Dep(true,    string(e.args[1])))
    end
    for arg in e.args
        if isa(arg, Expr)
            _extract_deps(arg, deps)
        end
    end
end

function extract_deps(e::Expr)
    deps = Dep[]
    _extract_deps(e, deps)
    return deps
end
