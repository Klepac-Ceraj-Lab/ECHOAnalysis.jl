# Metadata Analysis

## Getting Data

Now that we have the metadata in long form,
it's a bit easier to query it for the stuff we want.
I'm using the macros available from the [DataFramesMeta](https://github.com/JuliaData/DataFramesMeta.jl) package.

As an example, how many unique subjects do we have at least one sample for?

```@example 1
using ECHOAnalysis # hide
cd(joinpath(dirname(pathof(ECHOAnalysis)), "..")) # hide
```

```@example 1
using CSV
using DataFrames
using DataFramesMeta
using PrettyTables

allmeta = CSV.File("data/metadata/merged.csv") |> DataFrame

@linq allmeta |>
    select(:studyID) |>
    unique |> nrow
```

Or, how many subjects have two or more fecal samples?

```@example 1
sampleinfo = @linq allmeta |>
    where(:metadatum .== "CollectionRep") |>
    by(:studyID, nsamples = length(:studyID))

using Plots
using StatsPlots

histogram(sampleinfo[:nsamples], legend=false,
    title="Samples per Subject ID",
    xlabel="# of fecal samples", ylabel="# of subjects")

isdir("data/figures") || mkdir("data/figures")
savefig("data/figures/03-samples-per-subject.png") # hide
```

Wow - there are a couple of subjects that have a lot of samples.
Taking a look to see what's going on there:

```@example 1
# Which subjects are those?
highsamplers = @linq allmeta |>
    where(:metadatum .== "CollectionRep") |>
    by(:studyID, nsamples = length(:studyID)) |>
    where(:nsamples .>= 5) |>
    select(:studyID)
```

```@example 1
# filter on metagenomes (DOM = "Date of Metagenome")
@linq filter(r-> r[:studyID] in highsamplers[:studyID], allmeta) |>
    where(:metadatum .== "DOM") |>
    orderby(:studyID, :timepoint)

first(highsamplers, 15)
```

So a bunch of these are where
multiple samples were given for the same timepoint (eg `C0202_4F_1A` and `_2A`)
and/or both genotech and enthanol samples.


## Metagenomes

At this stage, what I care about are samples with metagenomes,
which are inidcated by the `DOM` metadatum.
To find all of the studyID/timepoint combos that have that have metagenomes:

```@example 1
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

filter!(mgxmeta) do row
    !any(ismissing.(row))
end

# show a random assortment of ~ 20 rows
mgxmeta[rand(nrow(mgxmeta)) .< 20 / nrow(mgxmeta), :] |> pretty_table
```

```@example 1
CSV.write("data/metadata/mgxmetadata.csv", mgxmeta);
```

Just checking that it reads in properly:

```@example 1
try
    mgxmeta = CSV.read("data/metadata/mgxmetadata.csv") |> DataFrame
    pretty_table(mgxmeta)
catch ex
    @error "Table didn't read properly" ex
end
```
