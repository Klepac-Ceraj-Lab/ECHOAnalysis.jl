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

const taxonlevels = BiobakeryUtils.taxonlevels

function widen2comm(df::DataFrame, features, samples; featurecol=:taxon, samplecol=:sample, abundcol=:abundance)
    @info "building feature dict"
    nf = length(features)
    featuredict = Dictionary(features, 1:nf)

    @info "building sample dict"
    ns = length(samples)
    sampledict = Dictionary(samples, 1:ns)

    @info "getting indices"
    idx = DataFrame(
        (sample  = sampledict[row[samplecol]],
        feature = featuredict[row[featurecol]])
        for row in eachrow(df))
    @info "filling matrix"
    ComMatrix(sparse(idx.feature, idx.sample, df[!,abundcol]), features, samples)
end

function taxonomic_profiles(taxprofile_path="/babbage/echo/profiles/taxonomic", taxlevel=:species;
                            filefilter=f->isstoolsample(basename(f)))
    filepaths = readdir(taxprofile_path, join=true)
    @info "Filtering files"
    filter!(filefilter, filepaths)
    df = DataFrame(sample=String[], taxon=String[], abundance=Float64[])

    features = Set(String[])
    samples  = Set(String[])

    @showprogress "Loading files" for file in filepaths
        sample = stoolsample(basename(file))
        push!(samples, sampleid(sample))
        tax = CSV.File(file, header=[:taxon, :taxid, :abundance, :additional_species],
                        skipto=5) |> DataFrame
        select!(tax, [:taxon, :abundance])
        taxfilter!(tax, taxlevel)
        union!(features, tax.taxon)

        tax.abundance ./= 100
        tax[!, :sample] .= sampleid(sample)
        append!(df, select(tax, [:sample, :taxon, :abundance]))
    end

    return df, features, samples
end

function functional_profiles(funcprofile_path="/babbage/echo/profiles/functional"; kind="ko_names_relab",
                            filefilter=f->isstoolsample(basename(f)), stratified=false)
    filepaths = readdir(funcprofile_path, join=true)
    @info "Filtering files"

    filter!(f-> occursin(kind, f) && filefilter(f), filepaths)
    df = DataFrame(sample=String[], func=String[], abundance=Float64[])

    features = Set(String[])
    samples  = Set(String[])

    @showprogress "Loading files" for file in filepaths
        sample = stoolsample(basename(file))
        push!(samples, sampleid(sample))
        func = CSV.File(file) |> DataFrame
        rename!(func, [:func, :abundance])

        stratified || filter!(row->!occursin("|g__", row.func) && !occursin("|unclassified", row.func), func)
        func[!,:sample] .= sampleid(sample)
        union!(features, func.func)
        append!(df, select(func, [:sample, :func, :abundance]))
    end

    return df, features, samples # widen2comm(df, featurecol=:func)
end


"""
    arrowconvert_taxonomic(filepaths, outpath)

"""
function arrowconvert_taxonomic(filepaths, outpath)  
    
end

function load_metaphlan_profile(file)
    sample = stoolsample(basename(file))
    sampledf = CSV.File(file, header=[:taxon, :taxid, :abundance, :additional_species],
                    datarow=5) |> DataFrame
    
    sampledf[!, :sample] .= sampleid(sample)
    sampledf.abundance ./= 100

    sampledf[!, :taxonlevel] .= :unidentified
    for rn in 1:nrow(sampledf)
        taxon, taxonlevel = last(parsetaxa(sampledf[rn, :taxon], throw=false))
        sampledf[rn, :taxon] = taxon
        sampledf[rn, :taxonlevel] = taxonlevel
    end

    select(sampledf, [:sample, :taxon, :taxonlevel, :taxid, :abundance])
end