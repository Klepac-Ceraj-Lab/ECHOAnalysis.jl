@info "Getting metadata"
include("script_helper.jl")

using SQLite
using ECHOAnalysis
using CSV

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "sqldb"
            help="path to SQLite database file"
        "--table", "-t"
            help="Table name in SQLite DB containing long form metadata"
            default="samples"
        "--metadata", "-m"
            help="""A list of metadata fields (columns) to include in wide table, separated by spaces.
                    Can be a path to a text file with one entry per line.
                    If this argument is not used, all metadata in the database will be included."""
            action=:append_arg
            nargs='+'
        "--parents", "-p"
            help="""A list of parent_table values containing metadata, separated by spaces.
                    Note: any fields passed to `metadata` that aren't from listed parents will be discarded.
                    Can be a path to a text file with one entry per line.
                    If this argument is not used, all parent tables in the database will be included."""
        action=:append_arg
        nargs='+'
        "--samples", "-s"
            help="""A list of samples separated by spaces, a path to a text file with one entry per line.
                    or a folder containing
                    If this argument is not used, all parent tables in the database will be included."""
            action=:append_arg
            nargs='+'
            required=true
        "--output", "-o"
            help="Path to save resulting table. Defaults to metadata.csv"
            default="metadata.csv"

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
    end

    return parse_args(s)
end


function main(args)

    if args["debug"]
        loglevel = Logging.Debug
    elseif args["verbose"]
        loglevel = Logging.Info
    elseif args["quiet"]
        loglevel = Logging.Error
    else
        loglevel = Logging.Warn
    end
    setup_logs!(loglevel, args["log"])

    @info args
    @info "loading db from $(args["sqldb"])"
    db = SQLite.DB(args["sqldb"])
    @info "getting long form metadata"
    longdf = getlongmetadata(db)

    parents = file_or_list_arg(args["parents"], longdf.parent_table, unique)
    metadata = file_or_list_arg(args["metadata"], longdf.metadatum, unique)

    samples = getsamples(args["samples"])

    @info "Widening"
    widedf = widemetadata(longdf, samples, parents=Set(parents), metadata=Set(metadata))
    @info "Saving wide metadata to $(args["output"])"
    CSV.write(expanduser(args["output"]), widedf)
end

# let args = parse_commandline()
#     main(args)
# end

# testing


let args = Dict{String,Any}("sqldb" => "/lovelace/echo/EchoSQL/databases/samples.sqlite",
                    "parents" => Array{Any,1}[],"output" => "~/Desktop/wide.csv",
                    "quiet" => false,"log" => nothing,"verbose" => true,
                    "metadata" => Array{Any,1}[],"table" => "samples","debug" => false,
                    "samples"=>"/lovelace/echo/sequencing/mgx/rawfastq/")
    main(args)
end
