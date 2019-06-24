function load_taxonomic_profiles(datatoml)
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

function load_metadata(datatoml, metadatakey="all"; samples::Union{Nothing,Vector{<:NamedTuple}}=nothing)
    long = CSV.read(datatoml["tables"]["metadata"][metadatakey]["path"])
    if !isnothing(samples)
        return getfocusmetadata(long, samples, focus=metadata_focus_headers)
    else
        return long
    end
end


const color1 = ColorBrewer.palette("Set1", 9)
const color2 = ColorBrewer.palette("Set2", 8)
const color3 = ColorBrewer.palette("Set3", 12)
const color4 = ColorBrewer.palette("Paired", 12)
const datatoml = parsefile("../../data/data.toml")
