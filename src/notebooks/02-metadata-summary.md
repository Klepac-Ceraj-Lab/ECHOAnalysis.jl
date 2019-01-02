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
    println(p)
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

parent_table = basename(metapaths[1])
println(parent_table)
```

```julia
fecaletoh = CSV.File(metapaths[1]) |> DataFrame

first(fecaletoh, 5)
```

Each table has it's own peculiarities,
but the end goal is to have a table where the columns are as follows:

| studyID        | timepoint           | date          | metadatum | value    | parent_table |
|----------------|---------------------|---------------|-----------|----------|--------------|
| Int - required | Int - if applicable | if applicable | required  | required | required     |
|                |                     |               |           |          |              |

For example:

```julia
for n in names(fecaletoh)
    s = string(n)
    if occursin("Fecalwethanol", s)
        s = replace(s, "Fecalwethanol"=>"")
    end
    s = replace(s, "#"=> "Num")
    rename!(fecaletoh, n=> Symbol(s))
end

using Dates

fecaletoh[:timestamp] = map(x-> DateTime(x, dateformat"m/d/y H:M:S"), fecaletoh[:timestamp])

fecaletoh[:date] = map(x-> DateTime(x, dateformat"m/d/y"),  fecaletoh[:CollectionDate])
deletecols!(fecaletoh, :CollectionDate)

fecaletoh = melt(fecaletoh, [:studyID, :date], variable_name=:metadatum)
fecaletoh[:parent_table] = "Fecal_with_Ethanol.csv"


# This table doesn't have "timepoint", so I'll add an empty column for now
# and eventually I'll get it from the date.

fecaletoh[:timepoint] = Vector{Union{Int, Missing}}(fill(missing, size(fecaletoh,1)))
```

It's mostly the same for Genotech ethanol samples.

```julia
parent_table = basename(metapaths[2])
println(parent_table)

fecal = CSV.File(metapaths[2]) |> DataFrame

fecal[:date] = map(x-> ismissing(x) ? missing : DateTime(x, dateformat"m/d/y"),  fecal[:collectionDate])
deletecols!(fecal, :collectionDate)

fecal = melt(fecal, [:studyID, :date], variable_name=:metadatum)
fecal[:parent_table] = parent_table

# This table doesn't have "timepoint", so I'll add an empty column for now
# and eventually I'll get it from the date.

fecal[:timepoint] = Vector{Union{Int, Missing}}(fill(missing, size(fecal,1)))
```

Now we can merge these tables together pretty easily

```julia
allmeta = vcat(fecal, fecaletoh)

@show first(allmeta, 3)
@show last(allmeta, 3)
```

Now for the other tables.

```julia
## Basic Family and Child Info
parent_table = basename(metapaths[3])
println(parent_table)

bfc = CSV.File(metapaths[3]) |> DataFrame

bfc = melt(bfc, :studyID, variable_name=:metadatum)

### In this case, the date and timepoint are irrelevant
bfc[:date] = Vector{Union{Dates.Date, Missing}}(fill(missing, size(bfc,1)))
bfc[:timepoint] = Vector{Union{Int, Missing}}(fill(missing, size(bfc,1)))
bfc[:parent_table] = parent_table

allmeta = vcat(allmeta, bfc)

## Delivery
parent_table = basename(metapaths[4])
println(parent_table)

delivery = CSV.File(metapaths[4]) |> DataFrame

delivery = melt(delivery, :studyID, variable_name=:metadatum)

### In this case, timepoint is irrelevant, but it would be nice to have the date.
### The closest we have is dueDate, but very few rows have it. Will leave missing for now
delivery[:date] = Vector{Union{Dates.Date, Missing}}(fill(missing, size(delivery,1)))
delivery[:timepoint] = Vector{Union{Int, Missing}}(fill(missing, size(delivery,1)))
delivery[:parent_table] = parent_table

allmeta = vcat(allmeta, delivery)

## Delivery
parent_table = basename(metapaths[5])
println(parent_table)

bfd = CSV.File(metapaths[5]) |> DataFrame

bfd = melt(bfd, :studyID, variable_name=:metadatum)

bfd[:date] = Vector{Union{Dates.Date, Missing}}(fill(missing, size(bfd,1)))
bfd[:timepoint] = Vector{Union{Int, Missing}}(fill(missing, size(bfd,1)))
bfd[:parent_table] = parent_table

allmeta = vcat(allmeta, bfd)

## Delivery
parent_table = basename(metapaths[6])
println(parent_table)

tpi = CSV.File(metapaths[6]) |> DataFrame
tpi[:date] = map(x-> ismissing(x) ? missing : DateTime(x, dateformat"m/d/y"),  tpi[:scanDate])
deletecols!(tpi, :scanDate)

tpi = melt(tpi, [:studyID, :date, :timepoint], variable_name=:metadatum)

tpi[:parent_table] = parent_table

allmeta = vcat(allmeta, tpi)

# show a random assortment of ~ 20 rows
@show allmeta[rand(size(allmeta, 1)) .< 20 / size(allmeta, 1), :]
```

```julia
CSV.write("data/metadata/merged.csv", allmeta)
```
