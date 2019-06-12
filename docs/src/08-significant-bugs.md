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
abtg = abundancetable(taxfilter(tax, :genus))
relativeabundance!(abt)
relativeabundance!(abtg)

kids_fecal = map(s-> occursin(r"^C", s[:sample]) && occursin("F", s[:sample]),
                    resolve_sampleID.(sitenames(abt)))

kidsp = view(abt, sites=kids_fecal)
kidgn = view(abtg, sites=kids_fecal)
```

## Define bugs

```@example sigbugs
genera = [
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
meta = getfocusmetadata(tables["metadata"]["merged_brain"]["path"],
    collect(samplenames(kidsp)))

meta[:sample] = [s.sample for s in resolve_sampleID.(samplenames(kidsp))]
meta = meta[[:subject,:timepoint,:sample, :correctedAgeDays]]

[collect(occurrences(view(kidgn, species=map(isequal("Bacteroides"),
    featurenames(kidgn)))))...]

for g in genera
    if g in featurenames(kidgn)
        @info g
        meta[Symbol(g)] = [collect(occurrences(view(kidgn, species=map(isequal("Bacteroides"),
            featurenames(kidgn)))))...]
    end
end

for s in species
    if s in featurenames(kidsp)
        @info s
        meta[Symbol(s)] = [collect(occurrences(view(kidsp, species=map(isequal(s),
            featurenames(kidsp)))))...]
    end
end


CSV.write("../../data/sigbugs.csv", meta)

metapa = deepcopy(meta)
for n in names(metapa[5:end])
    metapa[n] = Int.(present.(metapa[n]))
end

metapa
CSV.write("../../data/sigbugs_pa.csv", meta)

```
