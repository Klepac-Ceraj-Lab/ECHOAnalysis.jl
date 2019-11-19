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

function add_taxonomic_profiles(db::SQLite.DB, biobakery_path; replace=false, foldermatch=r"metaphlan2\/main")
    if "taxa" in DataFrame(SQLite.tables(db)).name
        !replace && error("Taxa already present in this database. Use `replace=true` to replace it")
        @warn "removing table taxa"
        SQLite.dropindex!(db, "samples_taxa")
        SQLite.drop!(db, "taxa")
    end

    taxlevels = BiobakeryUtils.taxlevels
    taxlevels = sort(collect(keys(taxlevels)), lt=(x,y)-> taxlevels[x] < taxlevels[y])
    for (root, dirs, files) in walkdir(biobakery_path)
        occursin(foldermatch, root) || continue
        filter!(f-> occursin(r"profile\.tsv", f), files)
        for file in files
            @info "Loading taxa for $file"
            sample = stoolsample(file)
            tax = CSV.read(joinpath(root, file), copycols=true)
            names!(tax, [:taxon, :abundance])

            for level in taxlevels
                filt = taxfilter(tax, level)

                total_abundance = sum(filt.abundance)
                filt.abundance ./= total_abundance

                filt[!, :kind] .= String(level)
                filt[!, :sample] .= sampleid(sample)
                SQLite.load!(filt, db, "taxa", ifnotexists=true)
            end
        end
    end
    @info "Creating sample index"
    SQLite.createindex!(db, "taxa", "samples_taxa", "sample", unique=false)
    @info "Done!"
end

function add_functional_profiles(db::SQLite.DB, biobakery_path;
        kind="genefamiles", stratified=false, replace=false,
        foldermatch=r"output\/humann2")
    if kind in DataFrame(SQLite.tables(db)).name
        !replace && error("$kind already present in this database. Use `replace=true` to replace it")
        @warn "removing table $kind"
        SQLite.dropindex!(db, "samples_$kind")
        SQLite.drop!(db, kind)
    end

    @info "Loading $kind functions, stratified = $stratified"
    for (root, dirs, files) in walkdir(biobakery_path)
        occursin(foldermatch, root) || continue
        filter!(f-> occursin(Regex(kind*".tsv"), f), files)

        for file in files
            @info "Loading functions for $file"
            sample = stoolsample(file)
            func = CSV.read(joinpath(root, file), copycols=true)
            names!(func, [:function, :abundance])
            func.stratified = map(row-> occursin(r"\|g__\w+\.s__\w+$", row[1]) ||
                                        occursin(r"\|unclassified$", row[1]),
                                        eachrow(func))

            stratified || filter!(row-> !row.stratified, func)

            func[!, :kind] .= kind
            func[!, :sample] .= sampleid(sample)
            SQLite.load!(func, db, kind, ifnotexists=true)
        end
    end

    @info "Creating sample index"
    SQLite.createindex!(db, kind, "samples_$kind", "sample", unique=false)
    @info "Done!"
end


function getlongmetadata(db::SQLite.DB, tablename="metadata")
    SQLite.Query(db, "SELECT * FROM '$tablename'") |> DataFrame
end


function makerowidx(df::DataFrame)
    Dict(r=>i for (i,r) in enumerate(df[!,1]))
end

function sqlprofile(db::SQLite.DB; tablename="taxa", kind="species", stratified=false, samplefilter=x->true)
    samples = stoolsample.(DataFrame(SQLite.Query(db, "SELECT DISTINCT sample FROM '$tablename'"))[!,1])
    filter!(samplefilter, samples)
    samples = Set(samples)

    @info "Creating table for $kind"
    cols = SQLite.columns(db, tablename)[!, :name]
    if stratified
        profile = SQLite.Query(db, "SELECT DISTINCT $(cols[1]) FROM '$tablename' WHERE kind='$kind' AND stratified=true") |> DataFrame
    else
        profile = SQLite.Query(db, "SELECT DISTINCT $(cols[1]) FROM '$tablename' WHERE kind='$kind'") |> DataFrame
    end

    ridx = makerowidx(profile)

    for s in samples
        @info "Loading $s"
        profile[!, Symbol(sampleid(s))] = Union{Float64,Missing}[missing for _ in eachrow(profile)]
        sdf = SQLite.Query(db, "SELECT $(cols[1]), abundance FROM '$tablename' WHERE kind='$kind' AND sample='$(sampleid(s))'") |> DataFrame

        profile[map(e-> ridx[e], sdf[!, Symbol(cols[1])]), Symbol(sampleid(s))] .= sdf[!,2]
    end

    @info "pruning empty rows"
    filter!(row-> any(!ismissing, row[2:end]), profile)
    replace!.(eachcol(profile[!,2:end]), Ref(missing=>0.))
    disallowmissing!(profile)
    return profile
end

sqlprofile(samplefilter, db::SQLite.DB; kwargs...) = sqlprofile(db; samplefilter=samplefilter, kwargs...)
