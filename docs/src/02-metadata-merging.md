# Working with Metadata

In the previous notebook, we generated separate metadata files
for different tables in the FilemakerPro database.
Now, we'll get these into a more usable format for analysis.

## Accessing TOML data in julia

Information about the locations of data are found in `data/data.toml`.
Parsing this file gives a set of nested key:value pairs.

```@example
cd(dirname(@__FILE__))
using ECHOAnalysis # hide
```

```@example 1
using Pkg.TOML: parsefile

files = parsefile("../../data/data.toml")
println.(keys(files))
files["description"]
```

The metadata tables are found under `["tables"]["metadata"]`

## Long form data

The `metascrubjl` script generates a table in long-form,
meaning each metadatum has its own row.


```@example 1
using CSV
using DataFrames
using PrettyTables

allmeta = CSV.File(files["tables"]["metadata"]["filemakerdb"]["path"]) |> DataFrame

# show ~10 random rows
allmeta[rand(nrow(allmeta)) .< 10 / nrow(allmeta), :] |> pretty_table
```

The script tries to find `timepoint` values for everything, and if it can't,
it assumes that the matadatum applies to all timepoints for that subject.
These are marked with `timepoint = 0`.
Let's look at which variables that applies to:

```@example 1
allmeta[allmeta[:timepoint] .== 0, :parent_table] |> unique;
```

Those all look reasonable!

## Sample Metadata

In addition to the FilemakerPro database,
we also have metdata info stored for each of the samples that are processed.
In this case, `timepoint` and `studyID` do not uniquely identify samples,
since we can have multiple samples per timepoint.
`sampleID`s should be unique though.

```@example 1
samples = CSV.File(files["tables"]["metadata"]["samples"]["path"]) |> DataFrame
rename!(samples, [:TimePoint=>:timepoint, :DOC=>:date, :SubjectID=>:studyID, :SampleID=>:sampleID])

samples = melt(samples, [:studyID, :timepoint, :sampleID], variable_name=:metadatum)
samples = filter(r-> !any(ismissing, [r[:studyID], r[:sampleID]]), samples)
disallowmissing!(samples)
samples[:parent_table] = "FecalProcessing";
```

We can only concatenate tables if they all have the same columns,
so I'll add a `sampleID` to all of the other observations.
The fecal sample `sampleID`s are build from the subjectID and timepoint,
so I'll do the same for other observations

```@example 1
allmeta[:sampleID] = map(r->
        "C" * lpad(string(r[:studyID]), 4, "0") * "_$(Int(floor(r[:timepoint])))M",
        eachrow(allmeta))

allmeta = vcat(allmeta, samples)
# reorder columns
allmeta = allmeta[[:sampleID, :studyID, :timepoint, :metadatum, :value, :parent_table]]
# remove rows with missing values
filter!(r-> !any(ismissing, [r[:studyID], r[:sampleID], r[:timepoint], r[:metadatum]]), allmeta)
disallowmissing!(allmeta)

# show a random assortment of ~ 10 rows
allmeta[rand(nrow(allmeta)) .< 10 / nrow(allmeta), :] |> pretty_table
```

A some of the `:value` fields have quotes or newlines in the field,
which screws up parsing later. For now I will just replace them.

```@example 1
for i in eachindex(allmeta[:value])
    v = allmeta[i, :value]
    ismissing(v) && continue
    if isa(v, AbstractString) && occursin(r"\n", v)
        v = replace(v, r"\n"=>"___")
        allmeta[i, :value] = v
    elseif isa(v, AbstractString) && occursin(r"\"", v)
        v = replace(v, r"\""=>"'")
        allmeta[i, :value] = v
    end
end


CSV.write("../../data/metadata/merged.csv", allmeta);
```
