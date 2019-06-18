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

ukids_abt = view(abt, sites=firstkids(samples))
ukids_dm = pairwise(BrayCurtis(), ukids_abt)
ukids_mds = fit(MDS, ukids_dm, distances=true)


allmeta = CSV.read("../../data/metadata/merged.csv")
filter(allmeta) do row
    row[:metadatum] == "motherSES"
end

focusmeta = getfocusmetadata("../../data/metadata/merged.csv",
    resolve_sampleID.(samplenames(ukids_abt)))
```


```@example glm
relativeabundance!(ukids_abt)

focusmeta[:sample] = samplenames(ukids_abt)

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

focusmeta[:breastfed] = breastfeeding.(eachrow(focusmeta))

@df focusmeta scatter(:correctedAgeDays ./ 365, :cogAssessment,
    ylabel="Composite score", xlabel="Age in years", legend=false)

@df focusmeta scatter(:white_matter_volume, :cogAssessment,
    ylabel="Composite score", xlabel="White Matter Volume", legend=false)

mylms = let df = DataFrame(bug=Symbol[], variable=String[], coef=Float64[], err=Float64[], t_value=Float64[], p_value=Float64[])
    for (i, sp) in enumerate(speciesnames(ukids_abt))
        sp = Symbol(sp)
        if prevalence(focusmeta[sp]) > 0.1
            mylm = @eval(lm(@formula($sp ~ cogAssessment + correctedAgeDays + motherSES + breastfed), focusmeta))
            tbl = coeftable(mylm)
            append!(df, DataFrame(
                bug      = sp,
                variable = tbl.rownms[2:end],
                coef     = coef(mylm)[2:end],
                err      = stderror(mylm)[2:end],
                t_value  = tbl.cols[3][2:end],
                p_value  = tbl.cols[4][2:end]
            ))
        end
    end
    df
end

using MultipleTesting
filter!(row-> !isnan(row[:p_value]) && row[:variable] != "correctedAgeDays", mylms)
mylms[:q_value] = adjust(mylms[:p_value], BenjaminiHochberg())
sort!(mylms, :q_value)
mylms[:q_value]
CSV.write("../../data/glms/julia-glms.csv", mylms)
```

```@example glm
oldermeta = filter(row-> !ismissing(row[:correctedAgeDays]) && row[:correctedAgeDays] / 365 > 1, focusmeta)

oldlms = let df = DataFrame(bug=Symbol[], variable=String[], coef=Float64[], err=Float64[], t_value=Float64[], p_value=Float64[])
    for (i, sp) in enumerate(speciesnames(ukids_abt))
        sp = Symbol(sp)
        if prevalence(oldermeta[sp]) > 0.1
            mylm = @eval(lm(@formula($sp ~ cogAssessment + correctedAgeDays + motherSES + breastfed), oldermeta))
            tbl = coeftable(mylm)
            append!(df, DataFrame(
                bug      = sp,
                variable = tbl.rownms[2:end],
                coef     = coef(mylm)[2:end],
                err      = stderror(mylm)[2:end],
                t_value  = tbl.cols[3][2:end],
                p_value  = tbl.cols[4][2:end]
            ))
        end
    end
    df
end

using MultipleTesting
filter!(row-> !isnan(row[:p_value]) && row[:variable] != "correctedAgeDays", oldlms)
oldlms[:q_value] = adjust(oldlms[:p_value], BenjaminiHochberg())
sort!(oldlms, :q_value)
oldlms[:q_value]
CSV.write("../../data/glms/julia-glms.csv", oldlms)

compl = focusmeta[completecases(focusmeta, [:correctedAgeDays, :cogAssessment, :motherSES, :breastfed]), :]

boxplot(compl[present.(compl[:Bacteroides_stercoris]), :motherSES], compl[present.(compl[:Bacteroides_stercoris]), :Bacteroides_stercoris])
scatter!(compl[:motherSES], compl[:Bacteroides_stercoris])
```


```@example glms
kids_spec = DataFrame(species=speciesnames(ukids_abt))

let sn = samplenames(ukids_abt)
    for i in eachindex(sn)
        kids_spec[Symbol(sn[i])] = occurrences(ukids_abt)[:, i]
    end
end

CSV.write("../../data/metadata/unique_kids_metadata.tsv",
            focusmeta[[:sample, :correctedAgeDays, :motherSES, :breastfed, :cogAssessment]],
            delim='\t')

!isdir("../../data/maaslin") && mkdir("../../data/maaslin")
CSV.write("../../data/maaslin/kids_species.tsv", kids_spec, delim='\t')

```
