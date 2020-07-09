using Distributed

MaxThreads = Sys.CPU_THREADS




# hire workers
if (length(workers())==1) && (workers()[1]==1)
    addprocs(MaxThreads)
elseif length(workers())<MaxThreads
    addprocs(MaxThreads-length(workers()))
end




# include the functions in each worker
@everywhere function LoadToThreadsLinux(i)
    printstyled("\nPlease wait while functions are loaded into each worker...\n", color=:green)
    tmp = remotecall(include, i, homedir() * "/Github/YujieGainRiskModel/main.jl")
    fetch(tmp)
end

@everywhere function LoadToThreadsWindows(i)
    printstyled("\nPlease wait while functions are loaded into each worker...\n", color=:green)
    tmp = remotecall(include, i,  homedir() * "\\Documents\\Github\\YujieGainRiskModel\\main.jl")
    fetch(tmp)
end

if VERSION>=v"1.0" && Sys.KERNEL==:Linux
    pmap(LoadToThreadsLinux, workers())
elseif VERSION>=v"1.0" && Sys.KERNEL==:NT
    pmap(LoadToThreadsWindows, workers())
else
    printstyled("\nThis version or system has not been tested...\n", color=:yellow)
end