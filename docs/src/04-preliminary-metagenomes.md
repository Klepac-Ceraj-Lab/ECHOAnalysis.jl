# Preliminary analysis of metagenomes

## Quality control

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

qc_files = filter(readdir(biobakery["kneaddata"]["path"])) do f
    occursin("read_counts", f)
end

qc_files = joinpath.(biobakery["kneaddata"]["path"], qc_files)
qc = CSV.File(qc_files[1]) |> DataFrame!
qc[:batch] = "batch001"

for f in qc_files[2:end]
    df = CSV.File(f) |> DataFrame!
    df[:batch] = match(r"(batch\d+)", f).captures[1]
    global qc = vcat(qc, df)
end

qc[:Sample] = map(qc[:Sample]) do s
    s = replace(s, "_kneaddata"=> "")
    resolve_sampleID(s)[:sample]
end
