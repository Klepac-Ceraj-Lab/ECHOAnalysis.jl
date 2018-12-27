using CSV
using DataFrames
using ArgParse
using Logging
using LoggingExtras


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
    "--in-folder", "-i"
        help="Folder containing fastq files"
        default="."
    "--delete"
        help="set to remove original files after concatenation"
        action=:store_true
    "--output-folder", "-o"
        help="destination folder for concatenated files"
        default="./"
    "--dry-run"
        help="print messages but do nothing"
        action=:store_true

    # Logging Stuff
    "--debug"
        help = "Show @debug level logging messages"
        action= :store_true
    "--verbose", "-v"
        help = "Show @info level logging messages"
        action= :store_true
    "--quiet", "-q"
        help = "Only show @error logging messages"
        action= :store_true
    "--log", "-l"
        help = "Write logs to this file. Default - no file"
        default = nothing
    end

    return parse_args(s)
end

function setup_logs!(loglevel, logpath; dryrun=false)
    glog = SimpleLogger(stderr, loglevel)
    if logpath === nothing || dryrun
        global_logger(glog)
    else
        global_logger(DemuxLogger(
            FileLogger(logpath, min_level=loglevel),
            glog, include_current_global=false))
    end
end

function catfiles(files::Vector{String}, outpath::String; keep=true, dryrun=false)
    @info "Concatenating:\n  $files\nto $outpath"
    if dryrun; @debug "Part of a dryrun, skipping"; return end

    open(outpath, "w") do outf
        for stream in open.(files)
            write(outf, stream)
            close(stream)
        end
    end
    (keep || dryrun) || rm.(files)
end

function main()
    args = parse_commandline()

    if args["debug"]
        loglevel = Logging.Debug
    elseif args["verbose"]
        loglevel = Logging.Info
    elseif args["quiet"]
        loglevel = Logging.Error
    else
        loglevel = Logging.Warn
    end
    setup_logs!(loglevel, args["log"], dryrun = args["dry-run"])

    if !args["dry-run"]
        isdir(args["output-folder"]) || mkdir(args["output-folder"])
    end

    k = !args["delete"]

    files = readdir(args["in-folder"])
    isdir("catfastq") || mkdir("catfastq")

    complete = String[]
    for f in files
        @debug f
        in(f, complete) && continue
        m = match(r"^([\w\-]+)_\w+?L00[0-9]_R(1|2)_001\.fastq\.gz", f)
        pid = m.captures[1]

        inds = findall(x-> occursin(Regex(pid), x), files)
        append!(complete, files[inds])

        pair1 = joinpath.(args["in-folder"], filter(x-> occursin(r"R1_001", x), files[inds]))
        pair2 = joinpath.(args["in-folder"], filter(x-> occursin(r"R2_001", x), files[inds]))

        catfiles(pair1,
            joinpath(args["output-folder"], "$pid.1.fastq.gz"),
            keep=k, dryrun=args["dry-run"])
        catfiles(pair2,
            joinpath(args["output-folder"], "$pid.2.fastq.gz"),
            keep=k, dryrun=args["dry-run"])
    end
end

main()
