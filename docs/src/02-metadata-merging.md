# Working with Metadata

In the previous notebook, we generated separate metadata files
for different tables in the FilemakerPro database.
Now, we'll get these into a more usable format for analysis.

## Accessing TOML data in julia

Information about the locations of data are found in `data/data.toml`.
Parsing this file gives a set of nested key:value pairs.

```@example metadata
cd(dirname(@__FILE__))
using ECHOAnalysis # hide
```

```@example metadata
using Pkg.TOML: parsefile

files = parsefile("../../data/data.toml")
println.(keys(files))
files["description"]
```

The metadata tables are found under `["tables"]["metadata"]`

## Long form data

The `metascrubjl` script generates a table in long-form,
meaning each metadatum has its own row.


```@example metadata
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

```@example metadata
allmeta[allmeta[:timepoint] .== 0, :parent_table] |> unique
```

Those all look reasonable!

## Sample Metadata

In addition to the FilemakerPro database,
we also have metdata info stored for each of the samples that are processed.
In this case, `timepoint` and `studyID` do not uniquely identify samples,
since we can have multiple samples per timepoint.
`sampleID`s should be unique though.

```@example metadata
samples = CSV.File(files["tables"]["metadata"]["samples"]["path"]) |> DataFrame
rename!(samples, [:TimePoint=>:timepoint, :DOC=>:date, :SubjectID=>:studyID, :SampleID=>:sampleID])

samples = melt(samples, [:studyID, :timepoint, :sampleID], variable_name=:metadatum)
samples = filter(r-> !any(ismissing, [r[:studyID], r[:sampleID]]), samples)
disallowmissing!(samples)
samples[:parent_table] = "FecalProcessing";
```

## Brain Data

Finally, we also have tables of brain volumes for many of our subjects.

```@example metadata
files["tables"]["brain"]["cortical"]["path"]
cortical = CSV.read(files["tables"]["brain"]["cortical"]["path"])
subcortical = CSV.read(files["tables"]["brain"]["subcortical"]["path"])


cortical[:studyID] = getfield.(parseletterid.(cortical[:SubjID]), :subject)
cortical[:timepoint] = getfield.(parseletterid.(cortical[:SubjID]), :timepoint)

# we only care about a subset of values for now
cortical = cortical[[:studyID, :timepoint, :LThickness, :RThickness,
                                           :LSurfArea,  :RSurfArea,
                                           :ICV]]

# convert to longform
cortical = melt(cortical, [:studyID, :timepoint], variable_name=:metadatum)
cortical[:parent_table] = "corticalVolumes"

# for the subcortex we mostly care about the total volume rather than individual values
subcortical[:subcortical_volume] = map(row-> sum(Vector(row[2:end-1])), eachrow(subcortical))
subcortical[:studyID] = getfield.(parseletterid.(subcortical[:SubjID]), :subject)
subcortical[:timepoint] = getfield.(parseletterid.(subcortical[:SubjID]), :timepoint)
subcortical = subcortical[[:studyID, :timepoint, :subcortical_volume]]

# convert to longform
subcortical = melt(subcortical, [:studyID, :timepoint], variable_name=:metadatum)
subcortical[:parent_table] = "subcorticalVolumes"
```

We can only concatenate tables if they all have the same columns,
so I'll add a `sampleID` to all of the other observations
to match what's in `samples`.
The fecal sample `sampleID`s are build from the subjectID and timepoint,
so I'll do the same for other observations

```@example metadata
allmeta[:sampleID] = map(r->
        "C" * lpad(string(r[:studyID]), 4, "0") * "_$(Int(floor(r[:timepoint])))M",
        eachrow(allmeta))

cortical[:sampleID] = map(r->
        "C" * lpad(string(r[:studyID]), 4, "0") * "_$(Int(floor(r[:timepoint])))M",
        eachrow(cortical))

subcortical[:sampleID] = map(r->
        "C" * lpad(string(r[:studyID]), 4, "0") * "_$(Int(floor(r[:timepoint])))M",
        eachrow(subcortical))


allmeta = vcat(allmeta, cortical, subcortical, samples)
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

```@example metadata
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
