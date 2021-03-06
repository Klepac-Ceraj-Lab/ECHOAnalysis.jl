using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))
push!(LOAD_PATH, @__DIR__)
@info "Scrubbing metadata"

using ArgParse
using Logging
using LoggingExtras

using CSV
using XLSX
using DataFrames
using Dates
using SQLite


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
            help = "for tabular text files, the delimeter used (generally ',' or '\t')"
            default = ","
        "--sheet"
            help = "for xlsx files, the name of the sheet that data is stored on"
            default = "Sheet1"
        "--dry-run"
            help = "Show logging, but take no action. Most useful with --verbose"
            action= :store_true
        "--output", "-o"
            help = "Output for scrubbed file (defaults to overwriting input)"
            arg_type = String
        "--force", "-f"
            help = "Allows overwriting of existing file. Reqired if no output path given"
            action=:store_true
        "--samples", "-s"
            help = "Path to samples to be included (optional). One line per sampleID."
            default = nothing
            arg_type = Union{Nothing, String}
        "--metadata", "-m"
            help = "Path to metada to be included (optional). One line per metadatum."
            default = nothing
            arg_type = Union{Nothing, String}
        "input"
            help = "Table to be scrubbed. By default, this will be overwritten if a CSV file"
            required = true
    end

    return parse_args(s)
end

function setup_logs!(loglevel, logpath; dryrun=false)
    glog = SimpleLogger(stderr, loglevel)
    if logpath === nothing || dryrun
        global_logger(glog)
    else
        logpath = abspath(expanduser(logpath))
        global_logger(
            TeeLogger(
                MinLevelLogger(FileLogger(logpath), loglevel),
                glog))
    end
end

function rowmissingfilter(r::DataFrameRow)
    if all(ismissing, values(r))
        @debug "All entries for row $(DataFrames.row(r)) are missing, removing"
        return false
    else
        return true
    end
end

function splitheader(h::AbstractString)
    h = replace(h, " "=>"_")
    return Tuple(String.(split(h, "::")))
end

function splitheader(h::Symbol)
    h = string(h)
    return splitheader(h)
end

function getsubtable(df, parent)
    indices = findall(x-> occursin(parent, replace(string(x), " "=>"_")), names(df))
    @debug "Getting subtable $parent at" indices
    table = df[:, indices]
    headers = map(x-> splitheader(x)[2], names(table))
    names!(table, Symbol.(headers), makeunique=true)

    return table
end


function elongate(df; idcol=:subject, tpcol=:timepoint)
    table = copy(df)
    if tpcol != :timepoint
        rename!(table, tpcol=>:timepoint)
    end

    table = melt(table, [idcol, tpcol], variable_name=:metadatum)
    filter!(r-> !ismissing(r[tpcol]), table)
    return table
end

function scrubdate!(df, colname)
    if eltype(df[!,colname]) <: AbstractString
        df.date = map(x-> ismissing(x) ? missing : DateTime(x, dateformat"m/d/y"),  df.date)
    elseif eltype(df[!, colname]) <: TimeType
        df.date = df[!, colname]
    end
    # delete original date column if it's not called `:date`
    colname != :date && select!(df, Not(colname))
end


# Use value types for special cases. See
#   https://docs.julialang.org/en/v1/manual/types/#%22Value-types%22-1
#   https://discourse.julialang.org/t/special-case-handling-alternatives-to-if-else/21608/4
struct ParentTable{T}
end

ParentTable(s::String) = ParentTable{Symbol(s)}()

# If no specific customprocess function exists, just return the table
customprocess!(table, ::ParentTable) = table

function customprocess!(table, ::ParentTable{:AAB})
    @warn "renaming subjectID"!
    rename!(table, [:subjectID=>:subject, :timePoint=>:timepoint])
end

function customprocess!(table, ::ParentTable{:Fecal_with_Ethanol})
    for n in names(table)
        s = string(n)

        # deal with headers like `Fecal_with_Ethanol::FecalwethanolShippedInitials`
        if occursin("Fecalwethanol", s)
            s = replace(s, "Fecalwethanol"=>"")
        end

        # deal with headers like `Fecal_with_Ethanol::Fecalwethanol#ofsamples`
        s = replace(s, "#"=> "Num")

        rename!(table, n=> Symbol(s))
    end
    rename!(table, :CollectionDate=>:date, :CollectionNum=>:timepoint)
    scrubdate!(table, :date)
    return table
end

function customprocess!(table, ::ParentTable{:FecalSampleCollection})
    rename!(table, :collectionDate=>:date, :collectionNum=>:timepoint)
    scrubdate!(table, :date)
    return table
end

function customprocess!(table, ::ParentTable{:TimepointInfo})
    scrubdate!(table, :scanDate)
    return table
end

function customprocess!(table, ::ParentTable{:GeneticOralSample})
    scrubdate!(table, :geneticCollectionDate)
    return table
end

function customprocess!(table, ::ParentTable{:GeneticOralSampleParent})
    scrubdate!(table, :CollectionDate)
    return table
end

function customprocess!(table,
            ::Union{ParentTable{:OralSampleCollection},ParentTable{:UrineSampleCollection}})
    scrubdate!(table, :collectionDate)
    return table
end



## Not sure this is the right thing to do - may nead to handle timepoint matching in separate script
# function customprocess!(table, ::ParentTable{:LeadHemoglobin})
#     rename!(table, :testNumber=>:timepoint)
#     return table
# end

function customprocess!(table, ::ParentTable{:Delivery})
    table.birthType = map(table.birthType) do t
        if ismissing(t)
            return missing
        elseif !occursin(r"([cC]esarean|[vV]aginal)", t)
            return t
        elseif occursin(r"[cC]esarean", t)
            return "Cesarean"
        elseif occursin(r"[vV]aginal", t)
            return "Vaginal"
        else
            @error "something went wrong" t
        end
    end
    return table
end

function loadtable(inputpath, delim, sheet)
    if endswith(inputpath, "xlsx")
        meta = DataFrame(
            XLSX.readtable(inputpath, sheet)...
            )
    else
        meta = DataFrame(
            CSV.File(inputpath, delim=delim)
            )
    end
    return meta
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
    setup_logs!(loglevel, args["log"], dryrun = args["dry-run"])

    inputpath = expanduser(args["input"]) |> abspath
    args["output"] === nothing ? outputpath = inputpath : outputpath = abspath(expanduser(args["output"]))
    (isfile(outputpath) && !args["force"]) && error("Output file already exists, used --force to overwrite")

    @info "Reading file from $inputpath"

    meta = loadtable(inputpath, args["delim"], args["sheet"])
    parents = map(n->splitheader(n)[1], names(meta)) |> unique

    @info "Processing parent tables" parents

    tables = DataFrame(subject=String[],
                        timepoint=Union{Int,Missing}[],
                        metadatum=String[],
                        value=Any[],
                        parent_table=String[])

    for p in parents
        subtable = getsubtable(meta, p)
        @info "Table Name: $p"

        # Remove columns where all values are missing
        for n in names(subtable)
            if all(ismissing, subtable[!, n])
                @info "All entries for subtable $p column $n are missing, removing"
                select!(subtable, Not(n))
            end
        end

        # Remove rows where all values are missing
        filter!(rowmissingfilter, subtable)

        # Rename columns to remove double :: and spaces
        for name in names(subtable)
            if occursin(" ", string(name))
                new_name = replace(string(name), " "=>"_")
                @info "Changing column $name to $new_name"
                rename!(subtable, name => new_name)
            end
        end

        customprocess!(subtable, ParentTable(p))
        rename!(subtable, :studyID=>:subject)

        if !any(n-> n == :timepoint, names(subtable))
            @warn "No timpoint column detected for $p, treating as all-timepoint variable"
            subtable[!, :timepoint] .= 0
        end
        nrow(subtable) < 2 && continue

        subtable = elongate(subtable, idcol=:subject, tpcol=:timepoint)
        subtable[!, :parent_table] .= p

        # One subject has a timepoint recorded as 2.5 for somet reason...
        if any(x-> x==2.5, subtable.timepoint)
            @warn "Removing timepoint 2.5 in table $p"
            filter!(r-> r.timepoint != 2.5, subtable)
            # convert column to Int
            subtable.timepoint = [t for t in subtable.timepoint]
        end

        tables = vcat(tables, subtable)
    end

    @info "Fixing weird strings"
    for (i, s) in enumerate(tables.value)
        isa(s, AbstractString) || continue
        s = replace(s, r"\n"=>"___")
        s = replace(s, r"\""=>"'")
        s = replace(s, r","=>";")
        tables[i, :value] = s
    end

    tables = unique(tables)

    if !isnothing(args["metadata"])
        keepers = Set(eachline(args["metadata"]))
        @info "Only keeping subset of metadata" keepers
        filter!(row-> string(row.metadatum) in keepers, tables)
    end

    @info "Writing scrubbed file to $outputpath"
    if !args["dry-run"]
        if endswith(outputpath, "sqlite")
            SQLite.DB(outputpath)
            SQLite.drop!(db, "metadata", ifexists=true)
            SQLite.load!(tables, db, "metadata")
            SQLite.dropindex!(db, "metadata_subjects_idx", ifexists=true)
            SQLite.createindex!(db, "metadata", "metadata_subjects_idx", "subject", unique=false)
            SQLite.createindex!(db, "metadata", "metadata_metadatum_idx", "metadatum", unique=false)
        else
            CSV.write(outputpath, tables, delim=args["delim"])
        end
    else
        @info "Just kidding! this is a dry run"
    end
end


let args = parse_commandline()
    main(args)
end
#
#
# ## Testing
#
# t = let args = Dict(
#             "input"    => "data/metadata/filemakerall.xlsx",
#             "output"   => "data/metadata/filemakerdb-restricted.csv",
#             "verbose"  => true,
#             "log"      => "data/metadata/filemaker-restricted.log",
#             "debug"    => false,
#             "quiet"    => false,
#             "delim"    => ",",
#             "dry-run"  =>false,
#             "sheet"    =>"Sheet1",
#             "metadata" => "data/metadata/keep_metadata.txt",
#             "force"    => true
#             )
#
#   main(args)
# end
