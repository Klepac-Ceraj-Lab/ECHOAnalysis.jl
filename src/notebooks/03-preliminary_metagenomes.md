# Preliminary analysis of metagenomes

## Quality control

```julia
using Pkg.TOML: parsefile
using CSV, DataFrames
using DataFramesMeta

biobakery = parsefile("data/data.toml")["tables"]["biobakery"]

kneadtables = joinpath.(biobakery["kneaddata"]["path"],
                      readdir(biobakery["kneaddata"]["path"]))

qc = vcat(map(t-> CSV.File(t, delim='\t') |> DataFrame, kneadtables))

names!(qc, map(n-> Symbol(replace(string(n), " "=> "_")), names(qc)))

println.(names(qc))

println.(qc[:Sample])

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
using Plots, StatPlots

let sids = Dict(y=>x for (x,y) in
    enumerate(unique(sort(@where(qc, :metadatum .== "raw"), :value)[:SampleID])))
    qc[:sid] = [sids[x] for x in qc[:SampleID]]
end



@df @where(qc, in.(:metadatum, Ref(["raw", "orphan", "final"]))
    ) groupedbar(:sid, :value, group=:metadatum, bar_position=:stack,
    xaxis="Samples", yaxis= "Count", legend=:topleft,
    title = "QC from Kneaddata"
    )
```



## Taxonomic profiles

We'll use some functions in the `Microbiome.jl` package

```julia
using Microbiome
using BiobakeryUtils
using MicrobiomePlots

mgxmeta = CSV.File("data/metadata/mgxmetadata.csv") |> DataFrame

@linq mgxmeta |>
    where(.!ismissing.(:SampleID), .!ismissing.(:value)) |>
    where(:metadatum .== "Mgx_batch") |>
    where((:value .== "Batch 1") .| (:value .== "Batch 2") .|
          (:value .== "Batch 3") .| (:value .== "Batch 4")) |>
    select(:SampleID) |>
    unique

samples = unique(skipmissing(mgxmeta[:SampleID]))
sids = unique(map(x-> String(split(x, "_")[1]), samples))



dfs = import_abundance_table.(joinpath.("data/biobakery/metaphlan2/",
        readdir("data/biobakery/metaphlan2/")))

samples = Set(Symbol[])
tax = DataFrame(SampleID=String[])

for (i, df) in enumerate(dfs)
    @info "batch $i"
    n = Set(names(df[2:end]))
    overlap = samples ∩ n
    deletecols!(df, [s for s in overlap])

    rename!(df, :col1=>:SampleID)
    global tax = join(tax, df, on=:SampleID, kind=:outer)

    global samples = samples ∪ n
end

for n in names(tax[2:end])
    tax[n] = Float64.(coalesce.(tax[n], 0.))
end



names!(tax, map(n-> replace(string(n), "_taxonomic_profile"=>"") |> Symbol, names(tax)))

dm = getdm(tax, BrayCurtis())
pco = pcoa(dm)

comm = abundancetable(tax)
featurenames(comm)

plot(pco, legend=false)
```
