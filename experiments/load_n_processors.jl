using Distributed

MaxThreads = Sys.CPU_THREADS




# hire workers
function LoadNProcessors(numb)
    # add workers
    if (length(workers())==1) && (workers()[1]==1) && (numb<=MaxThreads)
        addprocs(numb)
    elseif (length(workers())==1) && (workers()[1]==1) && (numb>MaxThreads)
        addprocs(MaxThreads)
    elseif length(workers())<numb && (numb<=MaxThreads)
        addprocs(numb-length(workers()))
    else
        addprocs(MaxThreads-length(workers()))
    end

    #include the functions in each worker
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
end
