function sampletable(rawfastq_path)
    ss = getsamples(rawfastq_path)
    df = DataFrame(ss)
    df.sample_type = sampletype.(ss)
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
            sample = resolve_sampleID(file)
            tax = CSV.read(joinpath(root, file), copycols=true)
            names!(tax, [:taxon, :abundance])

            for level in taxlevels
                filt = taxfilter(tax, level)

                total_abundance = sum(filt.abundance)
                filt.abundance ./= total_abundance

                filt[!, :kind] .= String(level)
                filt[!, :sample] .= sample.sample
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
