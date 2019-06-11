# Significant bugs

Looking at the literature and at the results from previous notebooks.

```@example sigbugs
cd(dirname(@__FILE__)) # hide
ENV["GKSwstype"] = "100" # hide

using ECHOAnalysis
using Pkg.TOML: parsefile
using DataFrames
using PrettyTables
using CSV
using Microbiome
using MultivariateStats
using StatsPlots
using MicrobiomePlots
using BiobakeryUtils
using ColorBrewer
using Clustering


# Setup

tables = parsefile("../../data/data.toml")["tables"]
figsdir = parsefile("../../data/data.toml")["figures"]["path"]
datafolder = tables["biobakery"]["path"]
metaphlan = tables["biobakery"]["metaphlan2"]
outdir = metaphlan["analysis_output"]
isdir(outdir) || mkdir(outdir)

tax = merge_tables(datafolder, metaphlan["root"], metaphlan["filter"],
    suffix="_profile.tsv")

# clean up sample names
names!(tax,
    map(n-> Symbol(
        resolve_sampleID(String(n))[:sample]),
        names(tax)
        )
    )

abt = abundancetable(taxfilter(tax))
pabt = abundancetable(taxfilter(tax, :phylum))
relativeabundance!(abt)
relativeabundance!(pabt)

kids = view(abt, sites=map(s-> occursin(r"^C", s[:sample]) && occursin("F", s[:sample]),
                    resolve_sampleID.(sitenames(abt))));
```

## Define bugs

```@example sigbugs
phyla = [
    "Bacillus",
    "Bacteroides",
    "Bifidobacterium",
    "Clostridiaceae",
    "Clostridiales",
    "Clostridium",
    "Enterrococcus",
    "Escherichia",
    "Lachnospiraceae",
    "Streptococcus",
    ]

species = [
    "Bacillus_subtilis",
    "Escherichia_coli",
    "Bifidobacterium_breve",
    "Bifidobacterium_longum",
    "Bifidobacterium_dentium",
    "Bifidobacterium_infantis",
    "Bifidobacterium_pseudocatenulatum"
]
```

## Sample metadata

```@example sigbugs

```
