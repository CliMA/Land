using Base.Test
using Autoreload

m1_before = """
module M
type MyType
  x::Int
end

f(var::MyType) = 1
end
"""

m1_after = """
module M
type MyType
  x::Int
end

f(var::MyType) = 2
end
"""

open(joinpath(dirname(@__FILE__), "M.jl"), "w") do f
    write(f, m1_before)
end

using Autoreload
arequire("M")
my_var = M.MyType(5)
M.f(my_var)

open(joinpath(dirname(@__FILE__), "M.jl"), "w") do f
    write(f, m1_after)
end

areload()
@test M.f(my_var) == 2

rm(joinpath(dirname(@__FILE__), "M.jl"))
