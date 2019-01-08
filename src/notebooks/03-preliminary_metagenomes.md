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
names!(tax, map(n-> replace(string(n), "-"=>"_") |> Symbol, names(tax)))

names!(tax, let x = []
        for name in String.(names(tax))
            try
                (pre, id, post) = match(r"([MC])(\d+)(_.+)", name).captures
                if length(id) == 3
                    id = "0" * id
                end
                push!(x, pre * id * post)
            catch e
                println(name)
                push!(x, name)
            end
        end
        Symbol.(x)
    end, makeunique=true)

taxfilter!(tax)

comm = abundancetable(tax)

dm = getdm(comm, BrayCurtis())
pco = pcoa(dm)

plot(pco, legend=false)
```

```julia

filter!(r-> all(!ismissing, [r[:studyID], r[:metadatum]]), mgxmeta)

samplemeta = let meta = Dict{String, Dict}()
    by(filter(r-> !ismissing(r[:SampleID]), mgxmeta), [:SampleID, :metadatum]) do df
        s = df[1,:SampleID]
        tp = df[1,:timepoint]
        id = df[1,:studyID]
        if !ismissing(s)
            if !haskey(meta, s)
                meta[s] = Dict{Any,Any}(:studyID => id, :timepoint=>tp)
            end

            m = Symbol(df[1,:metadatum])
            v = df[1,:value]
            nrow(df) > 1 && @warn ("too many rows for $s: $m") df

            if haskey(meta[s], m) && meta[s][m] != v
                error("duplicate entry for sample $s, metadatum $m, values: $(meta[s][m]) and $v")
            end

            meta[s][m] = v
        end
    end
    meta
end

by(filter(r-> ismissing(r[:SampleID]), mgxmeta), [:studyID, :timepoint, :metadatum]) do df
    id = df[1, :studyID]
    tp = df[1, :timepoint]

    entries = filter(k-> !ismissing(samplemeta[k][:studyID]) && samplemeta[k][:studyID] == id, keys(samplemeta))
    if !ismissing(tp)
        filter!(e-> !ismissing(samplemeta[e][:timepoint]) && samplemeta[e][:timepoint]== tp, entries)
    end

    m = df[1,:metadatum]

    if nrow(df) > 1
        u = unique(skipmissing(df[:value]))
        if length(u) == 0
            v = missing
        else
            length(u) > 1 && @warn ("too many rows for $id, $tp: $m") df
            v = u[1]
        end
    else
        v = df[1,:value]
    end

    for k in entries
        if haskey(samplemeta[k], m) && !ismissing(v)
            old = samplemeta[k][m]
            if !ismissing(old) && old != v
                if isa(old, Vector)
                    push!(samplemeta[k][m], v)
                else
                    samplemeta[k][m] = [old, v]
                end
                continue
            end
        end

        samplemeta[k][m] = v
    end
end

println.(sort(collect(keys(samplemeta))))
println.(samplenames(comm))

using ColorBrewer

c1 = ColorBrewer.palette("Set1", 9)

function getmeta(sampleid, metadatum, metadict)
    !haskey(metadict, sampleid) && return missing
    !haskey(metadict[sampleid], metadatum) && return missing
    return metadict[sampleid][metadatum]
end

function getcolors(vec, metadatum, metadict, colors)
    m = map(x-> getmeta(x, metadatum, metadict), vec)
    u = unique(m)
    length(u) > length(colors) && error("not enough colors")
    d = Dict(v => colors[i] for (i, v) in enumerate(u))
    return [d[x] for x in m]
end

function getcolors(vec, metadatum, metadict, colors::Dict)
    m = map(x-> getmeta(x, metadatum, metadict), vec)
    u = unique(m)
    length(u) > length(colors) && error("not enough colors")
    d = Dict(v => colors[i] for (i, v) in enumerate(u))
    return [d[x] for x in m]
end


csec = map(x-> getmeta(x, "birthType", samplemeta), samplenames(comm))
csec = map(x-> !ismissing(x) && occursin("Vaginal", x) ? "Vaginal" : x, csec)
csec = map(x-> !ismissing(x) && occursin("Cesarean", x) ? "Cesarean" : x, csec)
csec = map(x-> ismissing(x) ? "Missing" : x, csec)


plot(pco, group = csec)

plot(pco, group = map(x-> startswith(x, "M") ? "Mother" : "Child", samplenames(comm)))



age = map(x-> getmeta(x, "childCorrectedAgeMonths", samplemeta), samplenames(comm))
age = map(x-> ismissing(x) ? 150 / 12 : parse(Int, x) / 12, age)

plot(pco, zcolor=age)


CSV.write("problems.csv",
    DataFrame(sample_id=samplenames(comm),
             studyid = map(x-> getmeta(x, :studyID, samplemeta), samplenames(comm)),
             timepoint = map(x-> getmeta(x, :timepoint, samplemeta), samplenames(comm)),
             age = map(x-> getmeta(x, "childCorrectedAgeMonths", samplemeta), samplenames(comm))),
             )



```
