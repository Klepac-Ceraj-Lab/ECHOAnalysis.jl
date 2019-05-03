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
println.(keys(files))
files["description"]
```

The metadata tables are found under `["tables"]["metadata"]`


## Long form data

The `metascrubjl` script generates a table in long-form,
meaning each metadatum has its own row.


```julia
using CSV
using DataFrames

allmeta = CSV.File(files["tables"]["metadata"]["filemakerdb"]["path"]) |> DataFrame

@show allmeta[rand(nrow(allmeta)) .< 10 / nrow(allmeta), :]
```

The script tries to find `timepoint` values for everything, and if it can't,
it assumes that the matadatum applies to all timepoints for that subject.
These are marked with `timepoint = 0`.
Let's look at which variables that applies to:

```julia
by(allmeta[allmeta[:timepoint] .== 0, [:metadatum, :parent_table]], :metadatum) do df
    println("$(df[1,:parent_table]) : $(df[1,:metadatum])")
end
```

Those all look reasonable!

## Sample Metadata

In addition to the FilemakerPro database,
we also have metdata info stored for each of the samples that are processed.
In this case, `timepoint` and `studyID` do not uniquely identify samples,
since we can have multiple samples per timepoint.
`sampleID`s should be unique though.

```julia
samples = CSV.File(files["tables"]["metadata"]["samples"]["path"]) |> DataFrame
rename!(samples, [:TimePoint=>:timepoint, :DOC=>:date, :SubjectID=>:studyID, :SampleID=>:sampleID])

samples = melt(samples, [:studyID, :timepoint, :sampleID], variable_name=:metadatum)
filter!(r-> !any(ismissing, [r[:studyID], r[:sampleID]]), samples)
samples[:parent_table] = "FecalProcessing"
```

We can only concatenate tables if they all have the same columns,
so I'll add a `sampleID` to all of the other observations.
The fecal sample `sampleID`s are build from the subjectID and timepoint,
so I'll do the same for other observations

```julia
allmeta[:sampleID] = map(r->
        "C" * lpad(string(r[:studyID]), 4, "0") * "_$(Int(floor(r[:timepoint])))M",
        eachrow(allmeta))

allmeta = vcat(allmeta, samples)
# show a random assortment of ~ 20 rows
@show allmeta[rand(nrow(allmeta)) .< 10 / nrow(allmeta), :]

allmeta = allmeta[[:sampleID, :studyID, :timepoint, :metadatum, :value, :parent_table]]
filter!(r-> !any(ismissing, [r[:studyID], r[:sampleID], r[:timepoint], r[:metadatum]]), allmeta)
```

A some of the `:value` fields have quotes or newlines in the field,
which screws up parsing later. For now I will just replace them.

```julia
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


CSV.write("data/metadata/merged.csv", allmeta)
```


## Getting Data

Now that we have the metadata in this form,
it's a bit easier to query it for the stuff we want.
I'm using the macros available from the [DataFramesMeta](https://github.com/JuliaData/DataFramesMeta.jl) package.

As an example, how many unique subjects do we have at least one sample for?

```julia
using CSV
using DataFrames
using DataFramesMeta

allmeta = CSV.File("data/metadata/merged.csv") |> DataFrame

@linq allmeta |>
    select(:studyID) |>
    unique |> nrow
```

Or, how many subjects have two or more fecal samples?

```julia
sampleinfo = @linq allmeta |>
    where(:metadatum .== "CollectionRep") |>
    by(:studyID, nsamples = length(:studyID))

using Plots
using StatsPlots

histogram(sampleinfo[:nsamples], legend=false,
    title="Samples per Subject",
    xlabel="# of fecal samples", ylabel="# of subjects")
```

Wow - there are a couple of subjects that have a lot of samples.
Taking a look to see what's going on there:

```julia
# Which subjects are those?
highsamplers = @linq allmeta |>
    where(:metadatum .== "CollectionRep") |>
    by(:studyID, nsamples = length(:studyID)) |>
    where(:nsamples .>= 5) |>
    select(:studyID)

highsamplers = @linq filter(r-> r[:studyID] in highsamplers[:studyID], allmeta) |>
        where(:metadatum .== "DOM") |>
        orderby(:studyID, :timepoint)

first(highsamplers, 10)
```

So a bunch of these are where timepoints have multiple aliquots,
and/or both genotech and enthanol samples.


## Metagenomes

At this stage, what I care about are samples with metagenomes,
which are inidcated by the `DOM` metadatum.
To find all of the studyID/timepoint combos that have that have metagenomes:

```julia
mgxsamples = @linq allmeta |>
    where(:metadatum .== "DOM") |>
    select(:studyID, :timepoint) |>
    unique

sort!(mgxsamples, :studyID)

mgxmeta = let pairs = Set(zip(mgxsamples[:studyID], mgxsamples[:timepoint])), sids = Set(mgxsamples[:studyID])
    filter(row-> (row[:studyID] in sids && row[:timepoint] == 0) || # this captures metadata that's not linked to timepoints
                 (row[:studyID], row[:timepoint]) in pairs,
                 allmeta)
end

sort!(mgxmeta, [:studyID, :timepoint])
# show a random assortment of ~ 20 rows
@show mgxmeta[rand(nrow(mgxmeta)) .< 20 / nrow(mgxmeta), :]

isdir("data/biobakery/mgxmetadata/") || mkdir("data/biobakery/mgxmetadata/")
CSV.write("data/biobakery/mgxmetadata/mgxmetadata.csv", mgxmeta)
```

Just checking that it reads in properly:

```julia
mgxmeta = CSV.read("data/biobakery/mgxmetadata/mgxmetadata.csv") |> DataFrame
```

Not all of this metadata is useful for the first batch of samplese.
In addition, I want to subset the metadata into that that's relevant for
all timpoints, vs that that's timepoint-specific.

```julia
widedf = unstack(mgxmeta, [:studyID, :timepoint], :metadatum, :value)
CSV.write("data/biobakery/mgxmetadata/mgxmetadatawide.csv", widedf)

# timpoint 0 for all timepoint data
tp0 = filter(r-> r[:timepoint] == 0, widedf)
# some relevant info
tp0 = tp0[[:studyID, :timepoint, :birthType, :childGender, :childGestationalPeriodWeeks, :motherSES, :exclusivelyNursed]]
CSV.write("data/biobakery/mgxmetadata/mgxtp0.csv", tp0)

# What percent of subjects have each metadatum?
for n in names(tp0)
    println("Column $n: $(sum(ismissing, tp0[n]) / nrow(tp0))")
end

# how many subjects have all metadata?
sum(completecases(tp0))
```

```julia
# Everything that's NOT timpepoint 0
tps = filter(r-> r[:timepoint] != 0, widedf)
# relevant metadata
tps = tps[[:studyID, :timepoint, :correctedScanAge, :correctedAgeDays, :leadLevel,
            :childWeight, :childHeight, :childHeadCircumference]]

CSV.write("data/biobakery/mgxmetadata/mgxtps.csv", tps)
```
