const color1 = ColorBrewer.palette("Set1", 9)
const color2 = ColorBrewer.palette("Set2", 8)
const color3 = ColorBrewer.palette("Set3", 12)
const color4 = ColorBrewer.palette("Paired", 12)
const datatoml = parsefile("data/data.toml")

function load_taxonomic_profiles()
    bakery = datatoml["tables"]["biobakery"]

    tax = merge_tables(
            bakery["path"],
            bakery["metaphlan2"]["root"],
            bakery["metaphlan2"]["filter"],
            suffix="_profile.tsv")

    # clean up sample names
    names!(tax,
        map(n->
            Symbol(resolve_sampleID(String(n))[:sample]),
            names(tax)))

    return tax
end

function first_kids_tax_profiles(taxlevel=:species)
    tax = load_taxonomic_profiles()
    taxfilter!(tax, taxlevel)
    abt = abundancetable(tax)
    relativeabundance!(abt)
    kids = view(abt, sites=firstkids(samplenames(abt)))
    return kids
end

function load_functional_profiles(kind="genefamilies")
    bakery = datatoml["tables"]["biobakery"]

    tax = merge_tables(
            bakery["path"],
            bakery["humann2"]["root"],
            bakery["humann2"]["filter"],
            suffix="_genefamilies.tsv")

    # clean up sample names
    names!(tax,
        map(n->
            Symbol(resolve_sampleID(String(n))[:sample]),
            names(tax)))

    return tax
end

function load_metadata(datatoml=datatoml, metadatakey="all"; samples::Union{Nothing,Vector{<:NamedTuple}}=nothing)
    long = CSV.read(datatoml["tables"]["metadata"][metadatakey]["path"])
    if !isnothing(samples)
        wide = getfocusmetadata(long, samples, focus=metadata_focus_headers)
        return
    else
        return long
    end
end

function notebookpaths!(notebook)
    num = lpad(string(notebook), 2, "0")
    outpath = datatoml["notebooks"][num]["output"]
    figurepath = datatoml["notebooks"][num]["figures"]
    isdir(outpath) || mkpath(outpath)
    isdir(figurepath) || mkpath(figurepath)
    return outpath, figurepath
end
