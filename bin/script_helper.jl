using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))
push!(LOAD_PATH, @__DIR__)

using ArgParse
using Logging
using LoggingExtras

function setup_logs!(loglevel, logpath; dryrun=false)
    glog = SimpleLogger(stderr, loglevel)
    if logpath === nothing || dryrun
        global_logger(glog)
    else
        logpath = abspath(expanduser(logpath))
        global_logger(DemuxLogger(
            MinLevelLogger(FileLogger(logpath), loglevel),
            glog, include_current_global=false))
    end
end


function file_or_list_arg(arg, default, f=x->x)
    if length(arg) == 0
        return f(default)
    elseif length(arg) == 1 && isfile(arg[1])
        return readlines(arg[1])
    else
        return arg
    end
end
