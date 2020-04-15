function should_symbol_recurse(var)
    taboo = [Module, AbstractString, Dict, Array, Tuple, DataType, Function]
    for datatype in taboo
        if isa(var, datatype)
            return false
        end
    end
    return true
end

function module_rewrite(m::Module, name::Symbol, value)
    eval(m, :($name=$value))
end

function switch_mods(vars, mod1, mod2)
    types = extract_types(mod1)
    for var in vars
        var_type = typeof(var)
        if !isa(var_type,    DataType)
            continue
        end
        for (old_type, type_name) in types
            if is_same_type_base(var_type, old_type)
                mod2_raw_type = mod2.(type_name)
                mod2_type = mod2_raw_type{var_type.parameters...}
                if is_identical_type(var_type, mod2_type)
                    unsafe_alter_type!(var, mod2_type)
                else
                    warn("Couldn't alter $var")
                end
            end
        end
    end
end

macro type_strip()
    if VERSION >= v"0.3.0-prerelease+925"
        quote
        end
    else
        esc(quote 
            if options[:strip_types]
                e_i = eval_includes(e)
                if name == :Main
                    create_module(m_tmp, e)
                    e_t = strip_types(Main, e_i, Main.(m_tmp))
                else
                    create_module(m_tmp, Expr(:block)) 
                    info_debug("creating temporary module $m_tmp")
                    create_module(name, e, Main.(m_tmp))
                    info_debug("stripping types")
                    e_t = strip_types(m, e_i, Main.(m_tmp).(name))
                end
            else
                e_t = e
            end        
        end)
    end
end

function reload_module(name, e)
    local m_tmp
    info_debug("reloading module $name")
    while true
        m_tmp = symbol("_m_tmp_$(rand(1:100000))") #todo must be better way to do this
        if !(m_tmp in names(Main, true))
            break
        end
    end
    m = nothing
    if name in names(Main, true)
        m = Main.(name)
        if isa(m, Module)
            @type_strip
            # if options[:strip_types]
            #     e_i = eval_includes(e)
            #     if name == :Main
            #         create_module(m_tmp, e)
            #         e_t = strip_types(Main, e_i, Main.(m_tmp))
            #     else
            #         create_module(m_tmp, Expr(:block)) 
            #         info_debug("creating temporary module $m_tmp")
            #         create_module(name, e, Main.(m_tmp))
            #         info_debug("stripping types")
            #         e_t = strip_types(m, e_i, Main.(m_tmp).(name))
            #     end
            # else
            #     e_t = e
            # end
            e_t = e
            #module_rewrite(Main, name, m)
            if e_t != nothing
                info_debug("evaluating in context of original module")
                for arg in e_t.args
                    eval(m, arg)
                end
                return
            end
        end
        if name == :Main
            warn("Type defined in Main changed. Unavoidable error will occur")
        else
            warn("Recreating module $name because of changed types")
        end
    end
    info_debug("initializing module $name")
    create_module(name, e)
    if m != nothing
        if options[:smart_types]
            info_debug("switching types in $name")
            switch_types(m, Main.(name))
        end
    end
end

function eval_includes(e_block)
    e_list = []
    for e in e_block.args
        if isa(e, Expr) && 
            e.head == :call && 
            e.args[1] == :include 
            local f
            if isa(e.args[2], AbstractString)
                f = e.args[2]
            elseif isa(e.args[2], Expr)
                try
                    f = eval(Main, e.args[2]) # dangerous to eval in Main
                catch err
                    warn("include statement could not be processed:\n $err")
                    f = nothing
                end
            end
            if f != nothing
                f_p = parse_file(f)
                for arg in f_p.args
                    push!(e_list, arg)
                end
            else
                push!(e_list, e)
            end
            #push!(e_list, f_p)
        elseif isa(e, Expr) && e.head == :module
            f_p = eval_includes(e.args[3])
            push!(e_list, Expr(:module, e.args[1], e.args[2], f_p))
        else
            push!(e_list, e)
        end
    end
    block_wrap(e_list)
end


block_wrap(e_set::Array) = Expr(:block, e_set...)
block_wrap(e::Expr) = e

function get_type_name(e::Expr)
    if e.head == :type
        name = e.args[2]
    elseif e.head == :abstract
        name = e.args[1]
    else
        error("Name of type could not be extracted: \n$e")
    end
    while isa(name, Expr)
        name = name.args[1]
    end
    name
end

function strip_types(m, e_block, m_new)
    e_stripped = Any[]
    for e in e_block.args
        do_strip = false
        if isa(e, Expr) && e.head in [:type, :abstract]
            name = get_type_name(e) #e.args[2]
            if name in names(m, true)
                t = m.(name)
                if isa(t, DataType)
                    if is_identical_type(t, m_new.(name))
                        do_strip = true
                    else
                        info_debug("Type $t changed")
                        if DEBUG
                            module_rewrite(Main, :old_t, t)
                            module_rewrite(Main, :new_t, m_new.(name))
                        end
                        return nothing
                    end
                end
            end
        end
        if !do_strip
            push!(e_stripped, e)
        end
    end
    block_wrap(e_stripped)
end

function create_module(name, e, base=Main)
    info_debug("creating module")
    if name == :Main
        eval(Main, e)
    else
        eval(base, Expr(:module, true, name, e))
    end
end

function switch_types(m_old, m_new)
    symbols = collect_symbols(Main)
    switch_mods(symbols, m_old, m_new)
end

function extract_modules(e_block)
    modules = Dict{Symbol,Any}()
    in_main = Any[]
    info_debug("beginning module extraction")
    for e in e_block.args
        # isa(e, Expr) || continue
        if isa(e, Expr) && e.head == :module
            name = e.args[2]
            m_e = e.args[3]
            modules[name] = m_e
        else
            push!(in_main, e)
        end
    end
    modules[:Main] = block_wrap(in_main)
    info_debug("finishing module extraction")
    return modules
end

function collect_symbols(m)
    vars = []
    _collect_symbols(m, vars, 1)
    return vars
end

function try_getfield(m, name)
    local var
    try
        var = m.(name)
    catch err
        #warn("Symbol collection warning: \n$err")
        return (nothing, false)
    end
    return (var, true)
end

function _collect_symbols(m, vars, depth)
    # @show depth
    if isa(m, Module)
        name_list = names(m, true)
    else
        name_list = names(m)
    end
    for name in name_list
        var, has_var = try_getfield(m, name)
        has_var || continue
        if var in vars
            continue
        end
        push!(vars, var)
        if should_symbol_recurse(var)
            _collect_symbols(var, vars, depth+1)
        end
        if isa(var, Array) || isa(var, Tuple)
            for item in var
                if should_symbol_recurse(item)
                    _collect_symbols(item, vars, depth+1)
                end
                push!(vars, item)
            end
        end
        if isa(var, Dict) #todo handle keys properly
            for (key, value) in var
                for item in (key, value)
                    if should_symbol_recurse(item)
                        _collect_symbols(item, vars, depth+1)
                    end
                    push!(vars, item)
                end
            end
        end
    end
    vars
end

function unsafe_alter_type!(x, T::DataType)
    ptr = convert(Ptr{@compat UInt64}, pointer_from_objref(x))
    ptr_array = pointer_to_array(ptr, 1)
    ptr_array[1] = pointer_from_objref(T)
    x
end

function alter_type(x, T::DataType, var_name="")
    fields = names(T)
    old_fields = names(typeof(x))
    local x_new
    try
        x_new = T()
        for field in fields
            if field in old_fields
                x_new.(field) = x.(field)
            end
        end
    catch
        args = Any[]
        for field in fields
            if field in old_fields
                arg = x.(field)
                push!(args, arg)
            else
                warn("Type alteration of $(var_name) could not be performed automatically")
                return nothing
            end
        end
        x_new =    T(args...)
    end
    return x_new
end

function extract_types(mod::Module)
    types = Dict{DataType,Symbol}()
    for var_name in names(mod, true)
        var, has_var = try_getfield(mod, var_name)
        has_var || continue
        if typeof(var)==DataType
            types[var] = var_name
        end
    end
    types
end

function is_identical_type(t1, t2; kwargs...)
    t1 == t2
end

function is_identical_type(t1::TypeVar, t2::TypeVar; kwargs...)
    t1.name == t2.name || return false
    if isa(t1.ub, Tuple)
        isa(t2.ub, Tuple) || return false
        length(t1.ub) == length(t2.ub) || return false
        for i in 1:length(t1.ub)
            is_identical_type(t1.ub[i], t2.ub[i]) || return false
        end
    else
        is_identical_type(t1.ub, t2.ub) || return false
    end
    return true
end


function is_identical_type(t1::DataType, t2::DataType; deep_check::Bool=true)
    if !deep_check
        # This is hacky.
        return t1.name.name == t2.name.name
    end
    if t1.name.name == t2.name.name &&
        length(names(t1)) == length(names(t2)) && 
        all(names(t1).==names(t2))  
        #sizeof(t1)==sizeof(t2) && 
        # t1.parameters==t2.parameters && #verify this
        # all(fieldoffsets(t1).==fieldoffsets(t2)) &&
        # t1.types == t2.types
        if length(t1.parameters) != length(t2.parameters)
            return false
        end
        for i in 1:length(t1.parameters)
            if !is_identical_type(t1.parameters[i], t2.parameters[i], deep_check=false)
                return false
            end
        end
        for i in 1:length(t1.types)
            if !is_identical_type(t1.types[i], t2.types[i], deep_check=false)
                return false
            end
        end
        return true
    else
        false
    end
end


function is_same_type_base(t1::DataType, t2::DataType)
    if t1.name == t2.name
        true
    else
        false
    end
end
