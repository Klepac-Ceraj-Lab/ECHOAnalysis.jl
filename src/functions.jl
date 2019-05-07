"""
Parse a sample name and return its components.

Sample
"""
function resolve_sampleID(sid::Union{AbstractString,Symbol})
    sid = String(sid)
    m = match(r"^(([CM])(\d+)[_\-](\d)([FEO])[_\-](\d[A-Z]))(_S\d{1,2})?$", sid)


    if m == nothing
        return (sample=sid, subject=nothing, timepoint=nothing)
    else
        sample = String(m.captures[1])
        subject = parse(Int, m.captures[3])
        timepoint = parse(Int, m.captures[4])
        return (sample=sample, subject=subject, timepoint=timepoint)
    end
end

"""
Add batch info to `merged/` folders.

    Example:
        `batch001/metaplan2/merged/metaphlan2_taxonomic_profiles.tsv`
        to
        `batch001/metaplan2/merged/batch001_metaphlan2_taxonomic_profiles.tsv`
"""
function add_batch_info(folder)
    batch = ""
    for (root, dirs, files) in walkdir(folder)
        if occursin(r"batch\d+$", root)
            batch = match(r"(batch\d+)", root).captures[1]
            @warn batch
        end

        if occursin(r"merged$", root)
            for file in files
                if !occursin("batch", file)
                    oldpath = joinpath(root, file)
                    newpath = joinpath(root, "$(batch)_$file")
                    @info "moving $oldpath to $newpath"
                    mv(oldpath, newpath)
                end
            end
        elseif occursin(r"kneaddata$", root)
            if "kneaddata_read_counts.txt" in files
                oldpath = joinpath(root, "kneaddata_read_counts.txt")
                newpath = joinpath(root, "$(batch)_kneaddata_read_counts.txt")
                @info "moving $oldpath to $newpath"
                mv(oldpath,newpath)
            end
        end
    end
end

"""
Merge tables with a given subject based on their first column

Usage:
    merge_tables("metaplan", "profiles")
"""
function merge_tables(dir, suffix)
    isdir(dir) || throw(ArgumentError, "$dir is not a Directory")
    df = DataFrame(col1=[])

    for f in filter(x-> occursin(suffix, x), readdir(dir))
        t = CSV.File(joinpath(dir, f)) |> DataFrame!
        rename!(t, names(t)[1]=> :col1)
        df = join(df,t, on=:col1, kind=:outer)
    end
    return df
end
