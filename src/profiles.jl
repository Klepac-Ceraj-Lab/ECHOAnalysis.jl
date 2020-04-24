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

const taxlevels = BiobakeryUtils.taxlevels

function widen2comm(df::DataFrame; samplecol=:sample, featurecol=:taxon, abundcol=:abundance)
    @info "getting unique samples"
    samples = unique(df[!, samplecol])
    nsamples = length(samples)
    sampledict = HashDictionary(samples, 1:nsamples)

    @info "getting unique features"
    features = unique(df[!, featurecol])
    nfeatures = length(features)
    featuredict = HashDictionary(features, 1:nfeatures)

    mat = spzeros(nfeatures, nsamples)
    @info "filling matrix"
    for row in eachrow(df)
        mat[featuredict[row[featurecol]], sampledict[row[samplecol]]] = row[abundcol]
    end
    ComMatrix(mat, features, samples)
end

function taxonomic_profiles(taxprofile_path="/babbage/echo/profiles/taxonomic", taxlevel=:species;
                            filefilter=f->isstoolsample(basename(f)))
    filepaths = readdir(taxprofile_path, join=true)
    @info "Filtering files"
    filter!(filefilter, filepaths)
    df = DataFrame(sample=String[], taxon=String[], abundance=Float64[])

    @showprogress for file in filepaths
        sample = stoolsample(basename(file))
        tax = CSV.read(file, copycols=true)
        rename!(tax, [:taxon, :abundance])

        taxfilter!(tax, taxlevel)

        tax.abundance ./= 100
        tax[!, :sample] .= sampleid(sample)
        append!(df, select(tax, [:sample, :taxon, :abundance]))
    end

    return df # widen2comm(df)
end

function functional_profiles(funcprofile_path="/babbage/echo/profiles/functional"; kind="ko_names_relab",
                            filefilter=f->isstoolsample(basename(f)) && occursin(kind, f), stratified=false)
    filepaths = readdir(funcprofile_path, join=true)
    @info "Filtering files"
    filter!(filefilter, filepaths)
    df = DataFrame(sample=String[], func=String[], abundance=Float64[])

    @showprogress for file in filepaths
        sample = stoolsample(basename(file))
        func = CSV.File(file) |> DataFrame
        rename!(func, [:func, :abundance])

        stratified || filter!(row->!occursin("|g__", row.func) && !occursin("|unclassified", row.func), func)
        func[!, :sample] .= sampleid(sample)
        append!(df, select(func, [:sample, :func, :abundance]))
    end
    @info "Widening"

    return df # widen2comm(df, featurecol=:func)
end
