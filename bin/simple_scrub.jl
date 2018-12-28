using ArgParse
using Logging
using LoggingExtras
using CSV
using DataFrames


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        # Logging Stuff
        "--debug", "-d"
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

        # Stuff to do
        "--delim"
            default = ","
        "--dry-run"
            help = "Show logging, but take no action. Most useful with --verbose"
            action= :store_true
        "--output", "-o"
            help = "write to this file instead of overwriting original"
            arg_type = String
        "input"
            help = "Table to be scrubbed. By default, this will be overwritten"
            required = true
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

function rowmissingfilter(r::DataFrameRow)
    if all(ismissing, values(r))
        @info "All entries for row $(DataFrames.row(r)) are missing, removing"
        return false
    else
        return true
    end
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

    meta = CSV.File(args["input"], delim=args["delim"]) |> DataFrame

    # Remove columns where all values are missing
    for n in names(meta)
        if all(ismissing, meta[n])
            @info "All entries for column $n are missing, removing"
            delete!(meta, n)
        end
    end

    # Remove rows where all values are missing
    filter!(rowmissingfilter, meta)

    # Rename columns to remove double :: and spaces
    for name in names(meta)
        new_name = replace(string(name), " "=>"_")
        new_name = replace(new_name, "::"=>"__") |> Symbol

        @info "Changing column $name to $new_name"
        if !args["dry-run"]
            rename!(meta, name => new_name)
        else
            @info("Just kidding! this is a dry run")
         end
    end

    args["output"] === Nothing ? out = args["input"] : out = args["output"]
    out = abspath(expanduser(out))
    @info "Writing scrubbed file to $out"
    if !args["dry-run"]
        CSV.write(out, meta, delim=args["delim"])
    else
        @info "Just kidding! this is a dry run"
    end
end

main()
