# Preliminary analysis of metagenomes

## Quality control

```julia
using Revise
using Pkg.TOML: parsefile
using CSV, DataFrames
using DataFramesMeta

tables = parsefile("data/data.toml")["tables"]
biobakery = tables["biobakery"]

qc = CSV.File(biobakery["kneaddata"]["path"]) |> DataFrame

names!(qc, map(n-> Symbol(replace(string(n), " "=> "_")), names(qc)))

qc = @linq qc |>
  transform(raw = :raw_pair1 .+ :raw_pair2,
            trimmed = :trimmed_pair1 .+ :trimmed_pair2,
            orphan = :final_orphan1 .+ :final_orphan2,
            final = :final_pair1 .+ :final_pair2,
            )

qc = melt(qc, :Sample, variable_name=:metadatum)
rename!(qc, :Sample=>:SampleID)
qc[:metadatum] = string.(qc[:metadatum])
```

```julia
using Plots, StatsPlots

let sids = Dict(y=>x for (x,y) in
    enumerate(unique(sort(@where(qc, :metadatum .== "raw"), :value)[:SampleID])))
    qc[:sid] = [sids[x] for x in qc[:SampleID]]
end



@df @where(qc, in.(:metadatum, Ref(["raw", "orphan", "final"]))
    ) bar(:sid, :value, group=:metadatum,
    xaxis="Samples", yaxis= "Count", legend=:topleft,
    title = "QC from Kneaddata"
    )

isdir("data/figures") || mkdir("data/figures")
savefig("data/figures/03-knead-qc.svg")
```



## Taxonomic profiles

We'll use some functions in the `Microbiome.jl` package

```julia
using Microbiome
using BiobakeryUtils
using MicrobiomePlots

tp0 = CSV.File(tables["mgxmetadata"]["tp0"]["path"]) |> DataFrame
tps = CSV.File(tables["mgxmetadata"]["tps"]["path"]) |> DataFrame

tax = import_abundance_table.(tables["biobakery"]["metaphlan2"]["path"], delim=',')
names!(tax, map(x-> Symbol(replace(string(x), "_taxonomic_profile"=>"")), names(tax)))
rename!(tax, names(tax)[1]=>:taxon)

taxfilter!(tax, :species)

comm = abundancetable(tax)

dm = getdm(comm, BrayCurtis())
pco = pcoa(dm)

plot(pco, legend=false)
```

```julia

samples = samplenames(comm)
samples_match = map(s-> match(r"([CM])(\d+)-(\d)", s), samples)
samples_id = [x === nothing ? (0,0) : parse.(Int, (x.captures[2], x.captures[3])) for x in samples_match]

metadict = Dict(parse.(Int, (x.captures[2], x.captures[3]))=>
                Dict{Symbol, Any}(:subject_type=>string(x.captures[1])) for x in samples_match if x !== nothing
                )

for r in eachrow(tps)
    global metadict
    sid = (r[:studyID], Int(r[:timepoint]))
    !in(sid, keys(metadict)) && continue

    for n in names(r)
        metadict[sid][n] = r[n]
    end
end

for r in eachrow(tp0)
    global metadict
    sid = r[:studyID]
    ks = filter(k-> k[1] == sid, keys(metadict))
    for n in names(r[3:end])
        for k in ks
            metadict[k][n] = r[n]
        end
    end
end

using ColorBrewer

set1 = ColorBrewer.palette("Set1", 9)
set2 = ColorBrewer.palette("Set2", 8)

plot(pco, group=[s == (0,0) ? "ZymoControl" : metadict[s][:subject_type] for s in samples_id],
    legend=:topleft, color=[set1[4] set1[3] set1[9]], title="Taxonomic profiles - all samples")
savefig("data/figures/03-pcoa-mothers-children.svg")

kids = view(comm, sites = [i for i in eachindex(samples_id) if samples_id[i] != (0,0) && metadict[samples_id[i]][:subject_type] == "C"])
kids_id = [parse.(Int, (x.captures[2], x.captures[3])) for x in map(s-> match(r"([CM])(\d+)-(\d)", s),
    samplenames(kids))]

kids_dm = getdm(kids, BrayCurtis())
kids_pco = pcoa(kids_dm)

plot(pco, group=[:birthType in keys(metadict[s]) && !ismissing(metadict[s][:birthType]) ? metadict[s][:birthType] : "Missing" for s in kids_id],
    legend=:topleft, color=[set1[2] set1[9] set1[4]], title="Taxonomic profiles - birth type")
savefig("data/figures/03-pcoa-birthType.svg")

plot(pco, group=[:childGender in keys(metadict[s]) && !ismissing(metadict[s][:childGender]) ? metadict[s][:childGender] : "Missing" for s in kids_id],
    legend=:topleft, color=[set2[1] set2[2] set2[8]], title="Taxonomic profiles - Gender")
savefig("data/figures/03-pcoa-childGender.svg")


ages = view(comm, sites=[i for i in eachindex(samples_id) if samples_id[i] !== (0,0) &&
                                                            haskey(metadict[samples_id[i]], :correctedAgeDays) &&
                                                            !ismissing(metadict[samples_id[i]][:correctedAgeDays])])
ages_id = [parse.(Int, (x.captures[2], x.captures[3])) for x in map(s-> match(r"([CM])(\d+)-(\d)", s),
    samplenames(ages))]

ages_dm = getdm(ages, BrayCurtis())
ages_pco = pcoa(ages_dm)

plot(ages_pco, zcolor = log.([metadict[s][:correctedAgeDays] / 365 * 12 for s in ages_id]),
    title="Taxonomic profiles - Age in months (log)", color=:plasma, primary=false)
savefig("data/figures/03-pcoa-ages.svg")
```
