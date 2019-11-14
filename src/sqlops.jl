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

function taxonomic_profiles(db, biobakery_path)
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

function functional_profiles(db, biobakery_path)

end

function getlongmetadata(db, tablename="metadata")
    SQLite.Query(db, "SELECT * FROM '$tablename'") |> DataFrame
end

function sqlprofile(db, tablename="taxa")
    samples = DataFrame(SQLite.Query(db, "SELECT DISTINCT sample FROM '$tablename'"))[!,1] |> sort
    taxa = SQLite.Query(db, "SELECT DISTINCT taxon FROM '$tablename' WHERE kind='species'") |> DataFrame
    taxa |> SQLite.load!(db, "temp", temp=true, ifnotexists=true)

    for s in samples[1:1]
        sdf = SQLite.Query(db, "SELECT taxon abundance FROM '$tablename' WHERE kind='species' AND sample='$s'") |> DataFrame
    end

    SQLite.drop!(db, "temp")
end
