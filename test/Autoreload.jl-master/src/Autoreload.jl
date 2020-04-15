module Autoreload

using Compat

export arequire, areload,  aoptions_set, aoptions_get, smart_reload
export @aimport, @ausing

include("constants.jl")
include("smart_types.jl")
include("files.jl")
include("dependencies.jl")

DEBUG = false


info_debug(msg) = DEBUG ? info(msg) : nothing

function toggle_debug()
    global DEBUG
    DEBUG = !DEBUG
end

function aoptions_set(;kwargs...)
    global options
    for (key, value) in kwargs
        option_used = false
        for (i, (opt_key, opt_value)) in enumerate(options)
            if key==opt_key
                options[key] = value
                option_used = true
            end
        end
        if !option_used
            error("Option $key not recognized")
        end
    end
    options
end

aoptions_get() = copy(options)

function remove_file(filename)
    pop!(files, filename)
end

function arequire(filename=""; command= :on, depends_on=UTF8String[])
    if isempty(filename)
        return collect(keys(files))
    end
    original_filename = filename
    filename = find_file(standarize(filename), constants=options[:constants])
    filename == nothing && error("File $original_filename not found")
    if command in [:on, :on_depends]
        was_reloading = false
        if filename in keys(files)
            was_reloading = files[filename].should_reload
            remove_file(filename)
        end
        if command == :on
            should_reload = true
        else
            should_reload = was_reloading
        end
        files[filename] = AFile(should_reload, reload_mtime(filename), UTF8String[])
        parsed_file = parse_file(filename)
        auto_depends = extract_deps(parsed_file)
        auto_depends = [_.name for _ in auto_depends] #todo respect reload flag
        for i in 1:length(auto_depends)
            if !isabspath(auto_depends[i])
                auto_depends[i] = joinpath(dirname(filename), auto_depends[i])
            end
        end
        depends_on = vcat(auto_depends, depends_on)
        for d in depends_on
            d = standarize(d)
            if !haskey(files, d)
                d = find_file(d, base_file = filename)
                if d == nothing
                    error("Dependent file not found")
                end
                arequire(d, command = :on_depends)
            end
            push!(files[filename].deps, d)
        end
        if files[filename].should_reload
            include(find_file(filename))
        end
    elseif command == :off
        if haskey(files, filename)
            remove_file(filename)
        end
    else
        error("Command $command not recognized")
    end
    return
end

function smart_reload(file; kwargs...)
    global suppress_warnings
    info_debug("Initiating smart reload")
    original_file = file
    if !isabspath(file)
        file = find_file(standarize(file))
        if file == nothing
            error("File $(original_file) not found")
        end
    end
    suppress_warnings = true
    cd(dirname(file)) do
        info_debug("parsing $file")
        parsed = parse_file(file; kwargs...)
        info_debug("extracting modules")
        module_paths = extract_modules(parsed)
        info_debug("looping")
        for (module_name, e) in module_paths
            info_debug("reloading $module_name")
            reload_module(module_name, e)
            info_debug("done reloading $module_name")
        end
    end
    suppress_warnings = false
    return
end

function try_reload(file; kwargs...)
    smart_reload(file; kwargs...)
    # try
    #     smart_reload(file; kwargs...)
    # catch err
    #     if options[:verbose_level] == :warn
    #         warn("Could not autoreload $(file):\n$(err)")
    #     end
    # end
end

function is_equal_dict(d1::Dict, d2::Dict)
    for key in keys(d1)
        if d1[key]!=d2[key]
            return false
        end
    end
    return true
end

function areload(command= :force; kwargs...)
    if command == :use_state
        if options[:state] == :off
            return
        end
    elseif command == :off
        options[:state] = :off
        return
    elseif command == :on
        options[:state] = :on
    end
    dependencies = get_dependency_graph()
    file_order = topological_sort(dependencies)
    should_reload = [filename=>false for filename in file_order]
    marked_for_mtime_update = UTF8String[]
    for (i, file) in enumerate(file_order)
        file_time = files[file].mtime
        if reload_mtime(file) > file_time
            should_reload[file] = true
            # files[file].mtime = reload_mtime(file)
            push!(marked_for_mtime_update, file)
        end
    end
    old_reload = copy(should_reload)
    while true
        for (child, parents) in dependencies
            for parent in parents
                if should_reload[parent]
                    should_reload[child] = true
                end
            end
        end
        if is_equal_dict(old_reload, should_reload)
            break
        end
        old_reload = copy(should_reload)
    end

    for file in file_order
        if should_reload[file] && files[file].should_reload
            try_reload(file; constants=options[:constants], kwargs...)
        end
    end

    for file in marked_for_mtime_update
        files[file].mtime = reload_mtime(file)
    end

    return
end

function areload_hook()
    areload(:use_state)
end

if isdefined(Main, :IJulia)
    IJulia = Main.IJulia
    try
        IJulia.push_preexecute_hook(areload_hook)
    catch err
        warn("Could not add IJulia hooks:\n$(err)")
    end
end

include("macros.jl")

end
