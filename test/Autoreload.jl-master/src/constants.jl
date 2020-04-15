type AFile
    should_reload::Bool
    mtime::Float64
    deps::Vector{UTF8String}
end

type Dep
    should_reload::Bool
    name::UTF8String
end
suppress_warnings = false
const files = Dict{UTF8String,AFile}()
const options = @compat Dict{Symbol,Any}(:constants=>false, :strip_types=>true, :smart_types=>true, :verbose_level=>:warn, :state=>:on)
# verbose_level = :warn
# state = :on
