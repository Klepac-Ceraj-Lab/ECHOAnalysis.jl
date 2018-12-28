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



```julia
using CSV
using DataFrames


```
