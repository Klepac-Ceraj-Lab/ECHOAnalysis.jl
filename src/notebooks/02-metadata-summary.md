# Working with Metadata

In the previous notebook, we generated separate metadata files
for different tables in the FilemakerPro database.
Now, we'll get these into a more usable format for analysis.

## Accessing TOML data in julia

Information about the locations of data are found in `data/data.toml`.
Parsing this file gives a set of nested key:value pairs.

```julia
using Pkg.TOML: parsefile

files = parsefile("data/data.toml")
keys(files)
```

The metadata tables are found under `["tables"]["metadata"]`

```julia
println.(k for k in keys(files["tables"]["metadata"]))
```



```julia
metapaths = [files["tables"]["metadata"][k]["path"]
                for k in keys(files["tables"]["metadata"])]

for p in metapaths
    println(metapaths)
end
```

## Long form data

At the moment, the data has one subject or timepoint per row (wide-form).
To make things a bit easier to reshape and combine,
I'm going to convert it to long-form,
where each observation (metadatum) is on one row.

But we need a few pieces of info associated with each metadatum for them to
have meaning. In particular, we need the date of the observation,
the ID of the subject, and the timepoint number (if applicable).

```julia
using CSV
using DataFrames

bfc = CSV.File(metapaths[1]) |> DataFrame

first(bfc, 5)
```

Functions to sepparate two part column headers into their consituents.


```julia
function splitheader(h::AbstractString)
    return Tuple(String.(split(h, "__")))
end

function splitheader(h::Symbol)
    h = string(h)
    return splitheader(h)
end



bfc = melt(df, [:BasicFamilyAndChild__studyID])
rename!(bfc, :BasicFamilyAndChild__studyID=>:studyID)

function splitvariable!(df::AbstractDataFrame)
    df[:parent_table] = map(x-> splitheader(x)[1], df[:variable])
    df[:metadatum] = map(x-> splitheader(x)[2], df[:variable])
    delete!(df, :variable)
end

bfc[:parent_table] = map(x-> splitheader(x)[1], bfc[:variable])
bfc[:metadatum] = map(x-> splitheader(x)[2], bfc[:variable])

first(bfc, 5)
```

{{Descriptions of differences between tables}}


```julia
fecal = CSV.File(metapaths[2]) |> DataFrame

fecaletoh = fecal[map(n-> startswith(string(n), "Fecal_with_Ethanol"), names(fecal))]
filter!(x-> !ismissing(x[:Fecal_with_Ethanol__studyID]), fecaletoh)
rename!(fecaletoh, :Fecal_with_Ethanol__studyID=>:studyID)


fecaletoh = melt(fecaletoh, [:studyID, :Fecal_with_Ethanol__collectionDate])
splitvariable!(fecaletoh)




fecal = fecal[map(n-> !startswith(string(n), "Fecal_with_Ethanol"), names(fecal))]
filter!(x-> !ismissing(x[:FecalSampleCollection__studyID]), fecal)
rename!(fecal, :FecalSampleCollection__studyID=>:studyID)

fecal = melt(fecal, [:studyID, :FecalSampleCollection__collectionDate])
splitvariable!(fecal)



```
