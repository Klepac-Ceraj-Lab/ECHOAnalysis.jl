# Linear Models

```@example glm
cd(dirname(@__FILE__)) # hide
ENV["GKSwstype"] = "100" # hide

using ECHOAnalysis
using Pkg.TOML: parsefile
using DataFrames
using PrettyTables
using CSV
using Statistics
using Distances
using Microbiome
using MultivariateStats
using StatsPlots
using MicrobiomePlots
using BiobakeryUtils
using ColorBrewer
using Clustering

tables = parsefile("../../data/data.toml")["tables"]
figsdir = parsefile("../../data/data.toml")["figures"]["path"]
datafolder = tables["biobakery"]["path"]
metaphlan = tables["biobakery"]["metaphlan2"]
outdir = "../../data/glms"
isdir(outdir) || mkdir(outdir)

tax = merge_tables(datafolder, metaphlan["root"], metaphlan["filter"], suffix="_profile.tsv")

# clean up sample names
names!(tax,
    map(n-> Symbol(
        resolve_sampleID(String(n))[:sample]),
        names(tax)
        )
    )

taxfilter!(tax, :species)

abt = abundancetable(tax)
```

```@example glm
samples = resolve_sampleID.(samplenames(abt))
samples = samples[firstkids(samples)]

ukids_abt = view(abt, sites=firstkids(samples))
ukids_dm = pairwise(BrayCurtis(), ukids_abt)
ukids_mds = fit(MDS, ukids_dm, distances=true)


allmeta = CSV.read("../../data/metadata/merged.csv")
filter(allmeta) do row
    row[:metadatum] == "motherSES"
end

focusmeta = getfocusmetadata("../../data/metadata/merged.csv", samples)
```


```@example glm
relativeabundance!(ukids_abt)

focusmeta[:sample] = samplenames(ukids_abt)

kids_spec = DataFrame(species=speciesnames(ukids_abt))

let sn = samplenames(ukids_abt)
    for i in eachindex(sn)
        kids_spec[Symbol(sn[i])] = occurrences(ukids_abt)[:, i]
    end
end

CSV.write("../../data/metadata/unique_kids_metadata.tsv",
            focusmeta[[:sample, :correctedAgeDays, :motherSES, :birthType, :white_matter_volume, :grey_matter_volume, :csf_volume]],
            delim='\t')

!isdir("../../data/maaslin") && mkdir("../../data/maaslin")
CSV.write("../../data/maaslin/kids_species.tsv", kids_spec, delim='\t')
```

```@example glm
using RCall

R"""
 library(Maaslin2)
 fit_data <- Maaslin2("../../data/maaslin/kids_species.tsv", "../../data/metadata/unique_kids_metadata.tsv", "../../data/maaslin2_spec/")
 """
```

```@example glm
using GLM

sn = speciesnames(ukids_abt)
occ = occurrences(ukids_abt)
for i in eachindex(sn)
    focusmeta[Symbol(sn[i])] = asin.(sqrt.(collect(occ[i, :])))
end

let proj = projection(ukids_mds)
    for i in 1:size(proj, 2)
        focusmeta[Symbol("PCo$i")] = proj[:, i]
    end
end

focusmeta[:bayleysComposite] = map(row->
    mean([row[:languageComposite], row[:motorComposite]]),
    eachrow(focusmeta))

focusmeta[:cogAssessment] = map(eachrow(focusmeta)) do row
    cogs = row[[
                :mullen_EarlyLearningComposite,
                :fullScaleComposite,
                :FSIQ_Composite,
                :bayleysComposite
            ]]
    all(ismissing, cogs) && return missing
    return Float64(collect(skipmissing(cogs))[1])
end

focusmeta[:cogAssessment] = Union{Missing, Float64}[focusmeta[:cogAssessment]...]
describe(focusmeta[:cogAssessment])

@df focusmeta scatter(:correctedAgeDays ./ 365, :cogAssessment,
    ylabel="Composite score", xlabel="Age in years", legend=false)

@df focusmeta scatter(:ginisimpson, :cogAssessment,
    ylabel="Composite score", xlabel="Gini Simpson", legend=false)


@df focusmeta scatter(:white_matter_volume, :cogAssessment,
    ylabel="Composite score", xlabel="White Matter Volume", legend=false)


```


```@example glm



```
