# Omnibus Tests


```@example 1
# cd(dirname(@__FILE__)) # hide
using Pkg.TOML: parsefile
using CSV
using DataFrames
using Microbiome
using BiobakeryUtils
using ECHOAnalysis

tables = parsefile("data/data.toml")["tables"]

meta = CSV.read(tables["metadata"]["sample_metadata"]["path"])

sum(map(a-> !ismissing(a) && a < 366, meta[:correctedAgeDays]))

tax = merge_tables("data/biobakery/metaphlan/", "_profile.tsv")
# clean up sample names
names!(tax,
    map(n-> Symbol(
        resolve_sampleID(String(n))[:sample]),
        names(tax)
        )
    )

taxfilter!(tax, :species)

abt = abundancetable(tax)
dm = getdm(abt, BrayCurtis())

samples = samplenames(abt)
subject_type = [startswith(s, "M") ? "Mother" : "Child" for s in samples]

perm = permanova(dm, subject_type)
perm[:feature] = "species"
perm[:variable] = "Subject Type"

```


```@example 1
kids = view(abt, sites=map(s-> occursin(r"^C", s[:sample]) && occursin("F", s[:sample]),
                            resolve_sampleID.(sitenames(abt))))

kids_dm = getdm(kids, BrayCurtis())

p = permanova(kids_dm, string.(meta[:subject]))
p[:feature] = "species"
p[:variable] = "StudyID"
perm = vcat(perm, p)

unique_kids = let
    subjects= []
    unique = Bool[]
    for sample in sitenames(kids)
        s = resolve_sampleID(sample)
        if !in(s[:subject], subjects)
            push!(subjects, s[:subject])
            push!(unique,true)
        else
            push!(unique,false)
        end
    end
    unique
end



ukids = view(kids, sites=unique_kids)

ukids_dm = getdm(ukids, BrayCurtis())
ukids_pco = pcoa(ukids_dm)

ukidsmeta = meta[unique_kids, :]

p = permanova(ukids_dm, ukidsmeta[:correctedAgeDays])
p[:feature] = "species"
p[:variable] = "Age"
perm = vcat(perm, p)

p = permanova(ukids_dm, ukidsmeta[:birthType])
p[:feature] = "species"
p[:variable] = "birthType"
perm = vcat(perm, p)

p = permanova(ukids_dm, ukidsmeta[:breastfed])
p[:feature] = "species"
p[:variable] = "breastfed"
perm = vcat(perm, p)

p = permanova(ukids_dm, ukidsmeta[:formulafed])
p[:feature] = "species"
p[:variable] = "formulafed"
perm = vcat(perm, p)

p = permanova(ukids_dm, ukidsmeta[:childGender])
p[:feature] = "species"
p[:variable] = "childGender"
perm = vcat(perm, p)

p = permanova(ukids_dm, ukidsmeta[:motherSES])
p[:feature] = "species"
p[:variable] = "motherSES"
perm = vcat(perm, p)

youngkids = view(ukids, sites=map(a-> !ismissing(a) && a < 365*2, ukidsmeta[:correctedAgeDays]))
youngkids_dm = getdm(youngkids, BrayCurtis())
youngkidsmeta = ukidsmeta[map(a-> !ismissing(a) && a < 365*2, ukidsmeta[:correctedAgeDays]), :]

p = permanova(youngkids_dm, youngkidsmeta[:breastfed])
p[:feature] = "species"
p[:variable] = "young kids breastfed"
perm = vcat(perm, p)

p = permanova(youngkids_dm, youngkidsmeta[:formulafed])
p[:feature] = "species"
p[:variable] = "young kids formulafed"
perm = vcat(perm, p)




filter!(r-> !ismissing(r[Symbol("Pr(>F)")]), perm)
CSV.write("data/permanovas.csv", perm)
```



```@example 1
relativeabundance!(ukids)
kids_spec = DataFrame(species=speciesnames(ukids))

let sn = samplenames(ukids)
    for i in eachindex(sn)
        kids_spec[Symbol(sn[i])] = occurrences(ukids)[:, i]
    end
end

ukidsmeta[:sample] = samplenames(ukids)
ukidsmeta = ukidsmeta[[:sample, names(ukidsmeta[1:end-1])...]]
CSV.write("data/metadata/unique_kids_metadata.tsv",
            ukidsmeta[[:sample, :correctedAgeDays, :motherSES, :birthType, :white_matter_volume, :grey_matter_volume, :csf_volume]], delim='\t')

!isdir("data/maaslin") && mkdir("data/maaslin")
CSV.write("data/maaslin/kids_species.tsv", kids_spec, delim='\t')

#
# R"""
# library(Maaslin2)
# fit_data <- Maaslin2("data/maaslin/kids_species.tsv", "data/metadata/unique_kids_metadata.tsv", "data/maaslin2_spec/")
# """
```


```@example 1
hasbirthtype = .!ismissing.(ukidsmeta[:birthType])
hasbirthtype |> sum

bkidsmeta = ukidsmeta[hasbirthtype, :]
bkids = view(ukids, sites=hasbirthtype)
bkids_dm = getdm(bkids, BrayCurtis())

p = permanova(bkids_dm, bkidsmeta[:birthType])
p[:feature] = "species"
p[:variable] = "birthType"
perm = vcat(perm, p)



p = permanova(ukids_dm, string.(ukidsmeta[:breastfed]), 100_000)
```


```@example 1
newbatches = CSV.File("/Users/ksb/Desktop/batch89.csv") |> DataFrame

newbatches = filter(newbatches) do row
    row[:Mother_Child] == "C"
end

samples = resolve_sampleID.(newbatches[:SampleID])
append!(samples, resolve_sampleID.(readlines("samples.csv")))
append!(samples, resolve_sampleID.(samplenames(abt)))

subjects = [s.subject for s in samples]
timepoints = [s.timepoint for s in samples]

metadata = ["correctedAgeDays", "childGender", "birthType",
            "exclusivelyNursed", "exclusiveFormulaFed", "lengthExclusivelyNursedMonths",
            "amountFormulaPerFeed", "formulaTypicalType", "milkFeedingMethods",
            "typicalNumberOfEpressedMilkFeeds", "typicalNumberOfFeedsFromBreast",
            "noLongerFeedBreastmilkAge", "ageStartSolidFoodMonths"]

allmeta = CSV.File("data/metadata/merged.csv") |> DataFrame
focusmeta = getmetadata(allmeta, subjects, timepoints, metadata)

function convert2num2(s)
    ismissing(s) && return missing
    typeof(s) <: Number && return s
    if typeof(s) <: AbstractString
        s = occursin(".", s) ? parse(Float64, s) : parse(Int, s)
    end
    return s
end

focusmeta[:correctedAgeDays] = convert2num.(focusmeta[:correctedAgeDays])
focusmeta[:lengthExclusivelyNursedMonths] = convert2num.(focusmeta[:lengthExclusivelyNursedMonths])
focusmeta[:typicalNumberOfFeedsFromBreast] = convert2num.(focusmeta[:typicalNumberOfFeedsFromBreast])
focusmeta[:typicalNumberOfEpressedMilkFeeds] = convert2num.(focusmeta[:typicalNumberOfEpressedMilkFeeds])
focusmeta[:noLongerFeedBreastmilkAge] = convert2num.(focusmeta[:noLongerFeedBreastmilkAge])
focusmeta[:ageStartSolidFoodMonths] = convert2num.(focusmeta[:ageStartSolidFoodMonths])
focusmeta[:amountFormulaPerFeed] = convert2num.(focusmeta[:amountFormulaPerFeed])

filter!(focusmeta) do row
    !ismissing(row[:subject]) && !ismissing(row[:timepoint])
end





focusmeta[:breastfed] = breastfeeding.(eachrow(focusmeta))
focusmeta[:formulafed] = formulafeeding.(eachrow(focusmeta))

focusmeta |> CSV.write("/Users/ksb/Desktop/merged_metadata.csv")


```
