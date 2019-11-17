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

function taxonomic_profiles(db, biobakery_path; replace=false)
    if "taxa" in DataFrame(SQLite.tables(db)).name
        !replace && error("Taxa already present in this database. Use `replace=true` to replace it")
        SQLite.drop!(db, "taxa")
    end

    taxlevels = BiobakeryUtils.taxlevels
    taxlevels = sort(collect(keys(taxlevels)), lt=(x,y)-> taxlevels[x] < taxlevels[y])
    for (root, dirs, files) in walkdir(biobakery_path)
        occursin(r"metaphlan2\/main", root) || continue
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
                SQLite.load!(filt, db, "taxa")
            end
        end
    end
end

function functional_profiles(db, biobakery_path, kind="genefamiles")
    if kind in DataFrame(SQLite.tables(db)).name
        !replace && error("$kind already present in this database. Use `replace=true` to replace it")
        SQLite.drop!(db, "$kind")
    end


end

function getlongmetadata(db, tablename="metadata")
    SQLite.Query(db, "SELECT * FROM '$tablename'") |> DataFrame
end


function makerowidx(df::DataFrame)
    Dict(r=>i for (i,r) in enumerate(df[!,1]))
end

function sqlprofile(db; tablename="taxa", kind="species", samplefilter=x->true)
    samples = stoolsample.(DataFrame(SQLite.Query(db, "SELECT DISTINCT sample FROM '$tablename'"))[!,1])
    filter!(samplefilter, samples)

    taxa = SQLite.Query(db, "SELECT DISTINCT taxon FROM '$tablename' WHERE kind='$kind'") |> DataFrame
    ridx = makerowidx(taxa)

    for s in samples
        taxa[!, Symbol(s)] = Union{Float64,Missing}[missing for _ in eachrow(taxa)]
        sdf = SQLite.Query(db, "SELECT taxon, abundance FROM '$tablename' WHERE kind='species' AND sample='$s'") |> DataFrame

        taxa[map(e-> ridx[e], sdf.taxon), Symbol(s)] .= sdf[!,2]
    end

    replace!.(eachcol(taxa[!,2:end]), Ref(missing=>0.))
    disallowmissing!(taxa)
    return taxa
end

sqlprofile(samplefilter, db; kwargs...) = sqlprofile(db; samplefilter=samplefilter, kwargs...)
