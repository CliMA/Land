function find_in_path(name::AbstractString; constants::Bool=true)
    isabspath(name) && return name
    isfile(name) && return abspath(name)
    base = name
    if endswith(name, ".jl")
        base = name[1:end-3]
    else
        name = string(base,".jl")
        isfile(name) && return abspath(name)
    end
    for prefix in [Pkg.dir(); LOAD_PATH]
        path = joinpath(prefix, name)
        if constants
            file_name = name
        else
            file_name = string(name[1:end-3], "_reload", ".jl")
            if !isfile(prefix, base, "src", file_name)
                #warn("Reload version of $name not found")
                file_name = name
            end
        end
        isfile(path) && return abspath(path)
        path = joinpath(prefix, base, "src", file_name)
        isfile(path) && return abspath(path)
        path = joinpath(prefix, name, "src", file_name)
        isfile(path) && return abspath(path)
    end
    return nothing
end

function standarize(filename)
    if endswith(filename, ".jl")
        filename
    else
        string(filename, ".jl")
    end
end

function find_file(filename; base_file=nothing, constants=true)
    path = find_in_path(filename, constants=constants)
    if path == nothing && base_file!=nothing
        base_file = Base.find_in_path(base_file)
        path = joinpath(dirname(base_file), filename)
        if !isfile(path)
            path = nothing
        end
    end
    return path
end

function reload_mtime(filename)
    path = find_file(filename)
    m = mtime(path)
    if m == 0.0
        warn("Could not find edit time of $filename")
    end
    m
end

function parse_file(filename; kwargs...)
    path = find_file(filename; kwargs...)
    handle = open(path)
    source = string("begin\n", readall(handle), "\n end")
    close(handle)
    parsed = parse(source)
    return parsed
end
