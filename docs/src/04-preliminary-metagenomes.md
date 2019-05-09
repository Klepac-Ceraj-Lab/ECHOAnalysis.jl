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

```@example 1
tax = merge_tables("data/biobakery/metaphlan/", "_profile.tsv")
# clean up sample names
names!(tax,
    map(n-> Symbol(
        resolve_sampleID(String(n))[:sample]),
        names(tax)
        )
    )

using Microbiome
using MicrobiomePlots
using BiobakeryUtils


taxfilter!(tax)
abt = abundancetable(tax)

dm = getdm(tax, BrayCurtis())
pco = pcoa(dm)

plot(pco, legend=false, alpha=0.4)
```

```@example 1
c = [startswith(x, "C") ? :red : :blue for x in samplenames(abt)]

p = plot(pco, legend=false, marker=3,
    color=c, primary=false)
```


## Functions

```@docs
resolve_sampleID
```
