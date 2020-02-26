"""
    sampletable(rawfastqs::AbastractVector{<:AbstractString})

Get table with sample names and sample types
from a rawfastq folder.
"""
function sampletable(rawfastqs::AbstractVector{<:AbstractString})
    ss = stoolsample.(rawfastqs)
    df = DataFrame((sample     => sampleid(s),
                    subject    => subjectid(s),
                    timepoint  => timepoint(s),
                    sampletype => sampletype(s)
                    ) for s in ss)
    return df
end

"""
    add_taxonomic_profiles(db::SQLite.DB, biobakery_path; replace=false, foldermatch=r"metaphlan2/main")

Add taxonomic profiles to an SQLite database.
"""
function add_taxonomic_profiles(db::SQLite.DB, biobakery_path; replace=false, foldermatch=r"metaphlan2\/main")
    if "taxa" in SQLite.tables(db).name
        !replace && error("Taxa already present in this database. Use `replace=true` to replace it")
        @warn "removing table taxa"
        SQLite.dropindex!(db, "taxa_samples_idx", ifexists=true)
        SQLite.dropindex!(db, "taxa_taxon_idx", ifexists=true)
        SQLite.dropindex!(db, "taxa_batch_idx", ifexists=true)
        SQLite.drop!(db, "taxa")
    end

    taxlevels = BiobakeryUtils.taxlevels
    taxlevels = sort(collect(keys(taxlevels)), lt=(x,y)-> taxlevels[x] < taxlevels[y])
    filepaths = Tuple[]
    @info "Filtering files"
    for (root, dirs, files) in walkdir(biobakery_path)
        occursin(foldermatch, root) || continue
        batch = match(r"batch\d{3}", root)
        batch = isnothing(batch) ? "unknown" : String(batch.match)
        filter!(f-> occursin(r"profile\.tsv", f), files)
        append!(filepaths, [(batch, joinpath(root, f)) for f in files])
    end

    @info "Loading into database"
    @showprogress for (batch, file) in filepaths
        sample = stoolsample(basename(file))
        tax = CSV.read(file, copycols=true)
        names!(tax, [:taxon, :abundance])

        for level in taxlevels
            filt = taxfilter(tax, level)

            total_abundance = sum(filt.abundance)
            filt.abundance ./= total_abundance

            filt[!, :kind] .= String(level)
            filt[!, :sample] .= sampleid(sample)
            filt[!, :batch] .= batch
            SQLite.load!(filt, db, "taxa", ifnotexists=true)
        end
    end

    @info "Creating sample index"
    SQLite.createindex!(db, "taxa", "taxa_samples_idx", "sample", unique=false)
    SQLite.createindex!(db, "taxa", "taxa_taxon_idx", "taxon", unique=false)
    SQLite.createindex!(db, "taxa", "taxa_batch_idx", "batch", unique=false)
    @info "Done!"
end


"""
    add_functional_profiles(db::SQLite.DB, biobakery_path;
        kind="genefamiles_relab", stratified=false, replace=false,
        foldermatch=r"output/humann2", samples=:all)

Add functional profiles to an SQLite database.
Expects `kind` to come just before `.tsv` in filenames, eg
`C0001_1E_1A_genefamilies_relab.tsv` should have `kind="genefamilies_relab"`.
"""
function add_functional_profiles(db::SQLite.DB, biobakery_path;
            kind="genefamiles_relab", stratified=false, replace=false,
            foldermatch=r"output\/humann2", samples=:all)
    if kind in SQLite.tables(db).name
        !replace && error("$kind already present in this database. Use `replace=true` to replace it")
        @warn "removing table $kind"
        SQLite.dropindex!(db, "$(kind)_samples_idx", ifexists=true)
        SQLite.dropindex!(db, "$(kind)_function_idx", ifexists=true)
        SQLite.dropindex!(db, "$(kind)_batch_idx", ifexists=true)
        SQLite.drop!(db, kind)
        SQLite.drop!(db, "distinct_functions")
        SQLite.drop!(db, "distinct_samples")
    end

    @info "Loading $kind functions, stratified = $stratified"
    filepaths = Tuple[]
    for (root, dirs, files) in walkdir(biobakery_path)
        occursin(foldermatch, root) || continue
        batch = match(r"batch\d{3}", root)
        batch = isnothing(batch) ? "unknown" : String(batch.match)
        filter!(f-> occursin(Regex(kind*".tsv"), f), files)
        append!(filepaths, [(batch, joinpath(root, f)) for f in files])
    end

    @info "Loading functions for files"
    @showprogress for (batch, file) in filepaths
        sample = stoolsample(basename(file))

        (samples == :all || sample in samples) || continue
        func = CSV.read(file, copycols=true)
        names!(func, [:function, :abundance])
        func.stratified = map(row-> occursin(r"\|g__\w+\.s__\w+$", row[1]) ||
                                    occursin(r"\|unclassified$", row[1]),
                                    eachrow(func))

        stratified || filter!(row-> !row.stratified, func)

        func[!, :kind] .= kind
        func[!, :sample] .= sampleid(sample)
        func[!, :batch] .= batch
        SQLite.load!(func, db, kind, ifnotexists=true)
    end

    @info "Creating sample index"
    SQLite.createindex!(db, kind, "$(kind)_samples_idx", "sample", unique=false)
    @info "Creating feature index"
    SQLite.createindex!(db, kind, "$(kind)_function_idx", "function", unique=false)
    @info "Making distinct feature/sample tables"
    SQLite.Query(db, "SELECT DISTINCT function FROM $kind") |>
        SQLite.load!(db, "distinct_functions")
    SQLite.Query(db, "SELECT DISTINCT sample FROM $kind") |>
        SQLite.load!(db, "distinct_samples")

    @info "Done!"
end


function getlongmetadata(db::SQLite.DB, tablename="metadata")
    SQLite.Query(db, "SELECT * FROM '$tablename'") |> DataFrame
end


"""
    sqlprofile(db::SQLite.DB;
                    tablename="taxa", kind="species",
                    stratified=false, samplefilter=x->true,
                    prevalencefilter=(0.,1.)

Get a taxonomic or functional profile from SQLite database.
Returns dictionaries mapping sample-> column number
and feature-> row number,
along with a `SparseArray` of values.
"""
function sqlprofile(db::SQLite.DB;
                    tablename="taxa", kind="species",
                    stratified=false, samplefilter=x->true)
    @info "Starting"
    samples = let query = SQLite.Query(db, "SELECT DISTINCT sample FROM '$tablename'")
        [stoolsample(s[:sample]) for s in query]
    end

    filter!(samplefilter, samples)
    sort!(samples)
    sampleids = sampleid.(samples)
    sidx = dictionary(k=>i for (i,k) in enumerate(sampleids))

    @info "Creating table for $kind"
    cols = SQLite.columns(db, tablename)[!, :name]

    query = "SELECT DISTINCT $(cols[1]) FROM '$tablename' WHERE kind='$kind'" # AND sample IN ($(SQLite.esc_id(sampleid.(samples))))"
    stratified && (query *= " AND stratified=true")

    @info "Finding relevant features"
    if "distinct_functions" in SQLite.tables(db)[!,1]
        features = [x[1] for x in SQLite.Query(db, "SELECT DISTINCT $(cols[1]) FROM distinct_functions")]
    else
        features = [x[1] for x in SQLite.Query(db, query)]
    end
    fidx = dictionary(k=>i for (i,k) in enumerate(features))

    profile = spzeros(length(fidx), length(sidx))
    @info "Building profile"
    @showprogress 1 "Getting samples" for s in sampleids
        col = sidx[s]
        for r in SQLite.Query(db, "SELECT $(cols[1]), abundance FROM $tablename WHERE kind='$kind' AND sample='$s'")
            (feature, value) = (r[1], r[2])
            row = fidx[feature]
            profile[row, col] = value
        end
    end
    cm = ComMatrix(profile, features, sampleids)
    return view(cm, species=vec(sum(occurrences(cm), dims=2) .!= 0.))
end

sqlprofile(samplefilter, db::SQLite.DB; kwargs...) = sqlprofile(db; samplefilter=samplefilter, kwargs...)

function getallsamples(sqlite_path="/babbage/echo/sqldatabases/metadata.sqlite", table="allmetadata")
    db = SQLite.DB(sqlite_path)
    samples = SQLite.Query(db, "SELECT DISTINCT sample FROM $table") |> v-> [r.sample for r in v]
    filter!(s-> !occursin(r"^[CM]\d+_\d+M$", s), samples)
    return stoolsample.(samples)
end

function getmgxmetadata(sqlite_path="/babbage/echo/sqldatabases/metadata.sqlite", table="allmetadata"; samples=:all)
    db = SQLite.DB(sqlite_path)
    if samples == :all
        samples = SQLite.Query(db, "SELECT DISTINCT sample FROM $table") |> v-> [r.sample for r in v]
        filter!(s-> !occursin(r"^[CM]\d+_\d+M$", s), samples)
    end

    mdf = widemetadata(db, table, samples)
    filter!(row-> !ismissing(row.Mgx_batch), mdf)
    return mdf
end
