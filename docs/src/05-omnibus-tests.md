# Omnibus Tests


```julia
using Pkg.TOML: parsefile
using CSV
using DataFrames
using Microbiome
using BiobakeryUtils

tables = parsefile("data/data.toml")["tables"]

tp0 = CSV.File(tables["mgxmetadata"]["tp0"]["path"]) |> DataFrame
tps = CSV.File(tables["mgxmetadata"]["tps"]["path"]) |> DataFrame

tps[:timepoint] = Int.(tps[:timepoint])

tax = import_abundance_table.(tables["biobakery"]["metaphlan2"]["path"], delim=',')
samples = String.(names(tax)[2:end])
kids = filter(n-> startswith(n, "C"), samples)

kids_id = Set(map(x-> parse(Int,match(r"^C(\d+)_", x).captures[1]), kids))

kids_tp0 = filter(r-> r[:studyID] in kids_id, tp0)


kids_tps = filter(r-> r[:studyID] in kids_id, tps)

CSV.write("/Users/ksb/Desktop/kids_tp0.csv", kids_tp0)
CSV.write("/Users/ksb/Desktop/kids_tps.csv", kids_tps)
CSV.write("/Users/ksb/Desktop/kids_tp0_missing.csv", filter(r-> any(ismissing, r[2:end]), kids_tp0))

names!(tax, map(x-> Symbol(replace(string(x), "_taxonomic_profile"=>"")), names(tax)))
rename!(tax, names(tax)[1]=>:taxon)

taxfilter!(tax, :species)

tax = tax[[1, (sortperm(string.(names(tax[2:end]))) .+ 1)...]]
println.(names(tax))

let colnames = names(tax[2:end])
    subcol = Set([])
    keepers = Symbol[]
    for i in eachindex(colnames)
        m = match(r"(M|C)(\d+)-(\d)(\w)", string(colnames[i]))
        m === nothing && error("$(colnames[i]) didn't match regex")
    end
end

deletecols!(tax, Symbol("Zymo-Control"))

comm = abundancetable(tax)
dm = getdm(comm, BrayCurtis())

samples = samplenames(comm)
samples_match = map(s-> match(r"([CM])(\d+)-(\d)", s), samples)
samples_id = [x === nothing ? (0,0) : parse.(Int, (x.captures[2], x.captures[3])) for x in samples_match]

metadict = Dict(parse.(Int, (x.captures[2], x.captures[3]))=>
                Dict{Symbol, Any}(:subject_type=>string(x.captures[1])) for x in samples_match if x !== nothing
                )

for r in eachrow(tps)
    global metadict
    sid = (r[:studyID], Int(r[:timepoint]))
    !in(sid, keys(metadict)) && continue

    for n in names(r)
        metadict[sid][n] = r[n]
    end
end

for r in eachrow(tp0)
    global metadict
    sid = r[:studyID]
    ks = filter(k-> k[1] == sid, keys(metadict))
    for n in names(r[3:end])
        for k in ks
            metadict[k][n] = r[n]
        end
    end
end



perm = permanova(dm, [startswith(n, "C") ? "Child" : "Mom" for n in samplenames(comm)])
perm[:feature] = "species"
perm[:variable] = "subject type"

```


```julia
kids = view(comm, sites = [i for i in eachindex(samples_id) if samples_id[i] != (0,0) && metadict[samples_id[i]][:subject_type] == "C"])
kids_id = [parse.(Int, (x.captures[2], x.captures[3])) for x in map(s-> match(r"([CM])(\d+)-(\d)", s),
    samplenames(kids))]

kids_dm = getdm(kids, BrayCurtis())
p = permanova(kids_dm, [:birthType in keys(metadict[s]) &&
                    !ismissing(metadict[s][:birthType]) ?
                    metadict[s][:birthType] :
                    "Missing" for s in kids_id])

p[:feature] = "species"
p[:variable] = "birthType"
perm = vcat(perm, p)

p = permanova(kids_dm, [:exclusivelyNursed in keys(metadict[s]) &&
                    !ismissing(metadict[s][:exclusivelyNursed]) ?
                    metadict[s][:exclusivelyNursed] :
                    "Missing" for s in kids_id])

p[:feature] = "species"
p[:variable] = "exclusivelyNursed"
perm = vcat(perm, p)

p = permanova(kids_dm, [:childGender in keys(metadict[s]) &&
                    !ismissing(metadict[s][:childGender]) ?
                    metadict[s][:childGender] :
                    "Missing" for s in kids_id])

p[:feature] = "species"
p[:variable] = "childGender"
perm = vcat(perm, p)

ses = [:motherSES in keys(metadict[s]) &&
                    !ismissing(metadict[s][:motherSES]) ?
                    metadict[s][:motherSES] :
                    9999 for s in kids_id]

slvs = let slvs = String[]
    for s in ses
        if s == 9999
            push!(slvs, "unknown")
        elseif s < 6
            push!(slvs, "low")
        elseif s <= 7
            push!(slvs, "med")
        else
            push!(slvs, "high")
        end
    end
    slvs
end

newdm = DistanceMatrix(kids_dm[slvs .!= "unknown", slvs .!= "unknown"], BrayCurtis())
newslvs = slvs[slvs .!= "unknown"]

p = permanova(newdm, newslvs)


p[:feature] = "species"
p[:variable] = "motherSES"
perm = vcat(perm, p)


filter!(r-> !ismissing(r[Symbol("Pr(>F)")]), perm)
CSV.write("data/permanovas.csv", perm)


highres = CSV.File(tables["brain"]["highres"]["path"]) |> DataFrame
println.(names(highres))
highres[:studyID]

for r in eachrow(highres)
    startswith(r[:studyID], "m") && continue
    startswith(r[:studyID], "p") && continue

    id = (parse(Int,r[:studyID]), r[:timepoint])
    any(ismissing, id) && continue
    if haskey(metadict, id)
        if haskey(metadict[id], :correctedAgeDays)
            r[:correctedAgeDays] != metadict[id][:correctedAgeDays] &&
                @warn """
                Corrected age discrepancy for Subject $(id[1]), timepoint $(id[2])
                    Filemaker says: $(metadict[id][:correctedAgeDays])
                    HighRes table says: $(r[:correctedAgeDays])
                    Difference = $(metadict[id][:correctedAgeDays] - r[:correctedAgeDays])
                    """
        end

        for n in names(r[3:end])
            s = string(n)
            occursin(s, "Voxels") && continue
            s = replace(s, " (T)"=>"")
            s = replace(s, " "=>"_")
            s = Symbol(s)
            metadict[id][s] = r[n]
        end
    end
end


keys(metadict)
println.(names(highres))

highres[Symbol("White Matter Volume")]

wm = [:White_Matter_Volume in keys(metadict[s]) &&
                    !ismissing(metadict[s][:White_Matter_Volume]) ?
                    metadict[s][:White_Matter_Volume] : missing for s in kids_id]

gm = [:Gray_Matter_Volume in keys(metadict[s]) &&
                    !ismissing(metadict[s][:Gray_Matter_Volume]) ?
                    metadict[s][:Gray_Matter_Volume] : missing for s in kids_id]

wmgm = wm ./ gm
sum(!ismissing, wmgm)


newdm = DistanceMatrix(kids_dm[.!ismissing.(wmgm), .!ismissing.(wmgm)], BrayCurtis())

p = permanova(newdm, filter(!ismissing, wmgm))


p[:feature] = "species"
p[:variable] = "wm / gm"
perm = vcat(perm, p)


filter!(r-> !ismissing(r[Symbol("Pr(>F)")]), perm)
CSV.write("data/permanovas.csv", perm)

sum(z-> z>0, occurrences(kids[17:17, :]))

maaslin_df = DataFrame(
    subjectID = samplenames(kids),
    age = [metadict[s][:correctedAgeDays] for s in kids_id],
    white_gray_ratio = wmgm,
    birthType = [:birthType in keys(metadict[s]) && !ismissing(metadict[s][:birthType]) ?
        metadict[s][:birthType] : missing for s in kids_id],
    exclusivelyNursed = [:exclusivelyNursed in keys(metadict[s]) && !ismissing(metadict[s][:exclusivelyNursed]) ?
        metadict[s][:exclusivelyNursed] : missing for s in kids_id],
    childGender = [:childGender in keys(metadict[s]) && !ismissing(metadict[s][:childGender]) ?
        metadict[s][:childGender] : missing for s in kids_id]
        )

CSV.write("data/maaslin2_df.tsv", maaslin_df, delim='\t', missingstring="NA")


CSV.write("data/biobakery/metaphlan2/kids_spec.tsv",
    tax[[1, [i + 1 for i in eachindex(samples_id) if samples_id[i] != (0,0) &&
        metadict[samples_id[i]][:subject_type] == "C"]...]],
        delim='\t')

using RCall
using CSV
using DataFrames

maaslin_df = CSV.File("data/maaslin2_df.tsv", delim='\t') |> DataFrame

maaslin_df[103,:]
kids_spec = CSV.File("data/biobakery/metaphlan2/kids_spec.tsv", delim='\t') |> DataFrame
for n in names(kids_spec[2:end])
    s = sum(kids_spec[n])
    kids_spec[n] ./= s
end

names!(kids_spec, [Symbol(replace(string(n), '-'=>'_')) for n in names(kids_spec)])
CSV.write("data/biobakery/metaphlan2/kids_spec.tsv", kids_spec, delim='\t')

kids_spec_filt = filter(row-> sum(x-> x > 0, row[2:end]) / (ncol(kids_spec) - 1) >0.15, kids_spec)
CSV.write("data/biobakery/metaphlan2/kids_spec_filt.tsv", kids_spec_filt, delim='\t')

filter(r-> in(r[:subjectID], String.(names(kids_spec))), maaslin_df)


R"""
library(Maaslin2)
library(cplm)
fit_data <- Maaslin2("data/biobakery/metaphlan2/kids_spec.tsv", "data/maaslin2_df.tsv", "data/maaslin2_spec/")
"""


kids_spec

samples = filter(r-> ismissing(r[:white_gray_ratio]), maaslin_df)[:subjectID]

matches = map(s-> match(r"(C|M)(\d+)-(\d)", s), samples)
CSV.write("missing_highres.csv", DataFrame(subject=[parse(Int, m.captures[2]) for m in matches],
          timepoint = [parse(Int, m.captures[3]) for m in matches]))
