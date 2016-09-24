#
# Correctness Tests
# Many more should be included.

fatalerrors = length(ARGS) > 0 && ARGS[1] == "-f"
quiet = length(ARGS) > 0 && ARGS[1] == "-q"
anyerrors = false

using Base.Test
using SimpleGraphs

my_tests = ["test_cycles.jl", "test_transitiveclosure.jl"]

println("Running tests:")

for my_test in my_tests
    @show my_test
    try
        include(my_test)
        println("\t\033[1m\033[32mPASSED\033[0m: $(my_test)")
    catch e
        anyerrors = true
        println("\t\033[1m\033[31mFAILED\033[0m: $(my_test)")
        if fatalerrors
            rethrow(e)
        elseif !quiet
            showerror(STDOUT, e, backtrace())
            println()
        end
    end
end

if anyerrors
    throw("Tests failed")
end

#= I dont't know what it is for, so not included.
stdin = joinpath(dirname(@__FILE__), "stdin.sh")
ENV2 = copy(ENV)
ENV2["JULIA_HOME"] = JULIA_HOME
run(setenv(`bash $stdin`, ENV2))
=#
