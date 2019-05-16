# Preliminary analysis of metagenomes

All of the metagenomes were processed
using tools from the [bioBakery](https://bitbucket.org/biobakery/biobakery/wiki/Home).

```@example 1
using ECHOAnalysis
cd(joinpath(dirname(pathof(ECHOAnalysis)), "..")) # hide
```

```@example 1
using Pkg.TOML: parsefile
using CSV, DataFrames
using DataFramesMeta
using PrettyTables

tables = parsefile("data/data.toml")["tables"]
biobakery = tables["biobakery"]
```

## Quality control

First, I'll look at the QC results from `kneaddata`.

```@example 1
qc_files = filter(readdir(biobakery["kneaddata"]["path"])) do f
    occursin("read_counts", f)
end

# get all of the qc files
qc_files = joinpath.(biobakery["kneaddata"]["path"], qc_files)
# make a DataFrame of the first one
qc = CSV.read(qc_files[1])
qc[:batch] = "batch001"

# loop through the rest and concatenate to the first one
for f in qc_files[2:end]
    df = CSV.read(f)
    # get batch name from file
    df[:batch] = match(r"(batch\d+)", f).captures[1]
    global qc = vcat(qc, df)
end

pretty_table(qc)
```

To keep the formatting of sample IDs consistant across data types,
I'll use the [`resolve_sampleID`]@ref function.


```@example 1
qc[:Sample] = map(qc[:Sample]) do s
    s = replace(s, "_kneaddata"=> "")
    resolve_sampleID(s)[:sample]
end

pretty_table(qc)
```

```@example 1
# cleaning up the column names a bit
names!(qc,
    map(n-> Symbol(replace(String(n), " "=>"_")),
    names(qc)))
pretty_table(qc)
```

I don't really care about each mate pair individually,
so I'll sum them up

```@example 1
qc = @linq qc |>
  transform(raw = :raw_pair1 .+ :raw_pair2,
            trimmed = :trimmed_pair1 .+ :trimmed_pair2,
            orphan = :final_orphan1 .+ :final_orphan2,
            final = :final_pair1 .+ :final_pair2,
            )

```

```@example 1
using StatsPlots

sort!(qc, [:batch, :raw])

bar(x=qc[:Sample], hcat(qc[:raw], qc[:final]),
    xaxis="Samples", yaxis= "Count", legend=:topleft,
    title = "QC from Kneaddata", label=["Raw" "Final"],
    line=0)

savefig("data/figures/03-knead-qc.svg") # hide
```

These are a little more variable than I'd like.

```@example 1
using Statistics

qc_stats = by(qc, :batch) do df
                DataFrame(
                  mean=mean(df[:final]) / 1e6,
                  med=median(df[:final]) / 1e6,
                  max=maximum(df[:final]) / 1e6,
                  min=minimum(df[:final]) / 1e6,
                  )
end
```
```@example 1
CSV.write("data/biobakery/kneaddata/qc_stats.csv", qc_stats) # hide
```


## Taxonomic Profiles

Taxonomic profiles come from [MetaPhlAn2](https://bitbucket.org/biobakery/metaphlan2/src).
Each sample is run separately, and needs to be joined in a single table.
I'll use the function [`merge_tables`]@ref

```@example 2
using ECHOAnalysis
using DataFrames
using PrettyTables
using CSV
using Microbiome
using StatsPlots
using MicrobiomePlots
using BiobakeryUtils
using ColorBrewer

tax = merge_tables("data/biobakery/metaphlan/", "_profile.tsv")
# clean up sample names
names!(tax,
    map(n-> Symbol(
        resolve_sampleID(String(n))[:sample]),
        names(tax)
        )
    )
```

Some analysis of the fungi:

```@example 2
euk = filter(tax) do row
    occursin(r"^k__Eukaryota", row[1])
end

# remove columns that don't have any fungi
euk = euk[map(c->
    !(eltype(euk[c]) <: Number) || sum(euk[c]) > 0, names(euk))]

CSV.write("euk.csv", euk)
# get a df with only species
taxfilter!(euk)
CSV.write("euk_sp.csv", euk)
pretty_table(euk)
```

Those numbers are out of 100...
so really not much fungi at all,
at least according to metaplan.
There are some other methods to look more specifically at fungi,
which will have to wait for another time.

### PCoA Plots

For an initial overview,
let's look at the PCoA plots using BrayCurtis dissimilarity.

#### All Samples

```@example 2
spec = taxfilter(tax)
phyla = taxfilter(tax, :phylum)
spec |> pretty_table
```

```@example 2
abt = abundancetable(spec)
pabt = abundancetable(phyla)
relativeabundance!(abt)
relativeabundance!(pabt)

dm = getdm(spec, BrayCurtis())
pco = pcoa(dm)

plot(pco, legend=false, alpha=0.4)
```

```@example 2
color1 = ColorBrewer.palette("Set1", 9)
color2 = ColorBrewer.palette("Set2", 8)

c = [startswith(x, "C") ? color2[1] : color2[2] for x in samplenames(abt)]

p1 = plot(pco, marker=3, line=1, framestyle=1,
    color=c, primary=false)
scatter!([],[], color=color2[1], label="kids", legend=:topleft)
scatter!([],[], color=color2[2], label="moms", legend=:topleft)
title!("All samples taxonomic profiles")

savefig("data/figures/03-taxonomic-profiles-moms-kids.svg")
```

```@example 2
p2 = plot(pco, marker=3, line=1,
    zcolor=shannon(abt), primary = false, color=:plasma,
    title="All samples, shannon diversity")

savefig("data/figures/03-taxonomic-profiles-shannon.svg")
```

```@example 2
bacteroidetes = vec(Matrix(phyla[phyla[1] .== "Bacteroidetes", 2:end]))
firmicutes = vec(Matrix(phyla[phyla[1] .== "Firmicutes", 2:end]))

p3 = plot(pco, marker=3, line=1,
    zcolor=bacteroidetes, primary = false, color=:plasma,
    title="All samples, Bacteroidetes")

savefig("data/figures/03-taxonomic-profiles-bacteroidetes.svg")

p4 = plot(pco, marker=3, line=1,
    zcolor=firmicutes, primary = false, color=:plasma,
    title="All samples, Firmicutes")

savefig("data/figures/03-taxonomic-profiles-firmicutes.svg")

plot(p1, p2, p3, p4, marker = 2, markerstroke=0)
savefig("data/figures/03-taxonomic-profiles-grid.svg")
```

#### Kids

```@example 2
kids = view(abt, sites=startswith.(sitenames(abt), "C"))

kids_dm = getdm(kids, BrayCurtis())
kids_pco = pcoa(kids_dm)

plot(kids_pco, primary=false)
```

```@example 2
p5 = plot(kids_pco, marker=3, line=1,
    zcolor=shannon(kids), primary = false, color=:plasma,
    title="Kids, shannon diversity")

savefig("data/figures/03-taxonomic-profiles-kids-shannon.svg")

kids_bact = vec(collect(occurrences(view(pabt, species=occursin.("Bact", speciesnames(pabt))))))
kids_firm = vec(collect(occurrences(view(pabt, species=occursin.("Firm", speciesnames(pabt))))))
kids_act = vec(collect(occurrences(view(pabt, species=occursin.("Actino", speciesnames(pabt))))))
kids_proteo = vec(collect(occurrences(view(pabt, species=occursin.("Proteo", speciesnames(pabt))))))

plot(
    plot(kids_pco, marker=2, line=1,
        zcolor=kids_bact, primary = false, color=:plasma,
        title="Kids, Bacteroidetes"),
    plot(kids_pco, marker=2, line=1,
        zcolor=kids_firm, primary = false, color=:plasma,
        title="Kids, Firmicutes"),
    plot(kids_pco, marker=2, line=1,
        zcolor=kids_act, primary = false, color=:plasma,
        title="Kids, Actinobacteria"),
    plot(kids_pco, marker=2, line=1,
        zcolor=kids_proteo, primary = false, color=:plasma,
        title="Kids, Proteobacteria"),
    )
savefig("data/figures/03-taxonomic-profiles-kids-phyla.svg")
```
In order to decorate these PCoA plots with other useful information,
we need to return to the metadata.
I'll use the [`getmetadata`]@ref function


```@example 2
samples = resolve_sampleID.(samplenames(kids))

subjects = [s.subject for s in samples]
timepoints = [s.timepoint for s in samples]
metadata = ["correctedAgeDays","childGender","APOE","birthType","exclusivelyNursed","exclusiveFormulaFed","lengthExclusivelyNursedMonths","noLongerFeedBreastmilkAge","ageStartSolidFoodMonths","motherSES","childHeight","childWeight",]

df = getmetadata(allmeta, subjects, timepoints, metadata)
df[:motherSES] = map(x-> x == "9999" ? missing : parse(Int, x), df[:motherSES])

df |> CSV.write("focus2.csv") # hide
```

##### Birth type

```@example 2
plot(kids_pco, marker=3, line=1,
    color=metacolor(df[:birthType], color2[4:5], missing_color=color2[end]),
    title="Kids, BirthType", primary=false)
scatter!([],[], color=color2[4], label=unique(df[:birthType])[1])
scatter!([],[], color=color2[5], label=unique(df[:birthType])[2])
scatter!([],[], color=color2[end], label="missing")

savefig("data/figures/03-taxonomic-profiles-kids-birth.svg")
```

##### Breastfeeding

Information about braestfeeding is spread across 2 different parent tables.
`BreastfeedingDone` indicates that the child is no longer breastfeeding,
and has a lot of information about formula use, solid food etc,
`BreastfeedingStill` is for kids that are still breastfeeding,
and has substantially less information.

```
```

To make it a bit easier to get a handle on,
I created a wide-version table with only these parent tables,
sorted by subjectID.
The `BreastfeedingStill` values are not timepoint associated,
while the `BreastfeedingDone` are.


```@example 2
bf = filter(allmeta) do row
    row[:parent_table] == "BreastfeedingStill" || row[:parent_table] == "BreastfeedingDone"
end

bf = unstack(bf, [:studyID, :timepoint], :metadatum, :value)
bf |> CSV.write("data/metadata/breastfeeding.csv")
```

## Functions

```@docs
resolve_sampleID
merge_tables
getmetadata
```
