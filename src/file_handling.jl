"""
    add_batch_info(folder::AbstractString)

Add batch info to `merged/` folders.

Example:
    `batch001/metaplan2/merged/metaphlan2_taxonomic_profiles.tsv`
    to
    `batch001/metaplan2/merged/batch001_metaphlan2_taxonomic_profiles.tsv`
"""
function add_batch_info(folder::AbstractString)
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
    merge_tables(folder, dataroot, fltr; suffix=fltr)

Recurssively walk through a `folder` and merge tables inside folders named
`dataroot` with a given filter `fltr`. The `fltr` is assumed to be the
end of the filename and removed, but a different `suffix` may be provided.

Usage:
    merge_tables("biobakery", "metaphlan2", "_profile.tsv")
"""
function merge_tables(folder, dataroot, fltr; suffix=fltr)
    isdir(folder) || throw(ArgumentError("$folder is not a Directory"))
    df = DataFrame(col1=[])

    for (root, dirs, files) in walkdir(folder)
        occursin(dataroot, root) || continue

        for f in filter(x-> occursin(fltr, x), files)
            t = CSV.File(joinpath(root, f)) |> DataFrame!
            fname = Symbol(replace(f, suffix=>""))
            rename!(t, names(t)[1]=> :col1, names(t)[2]=>fname)
            df = join(df,t, on=:col1, kind=:outer)
        end
    end
    for n in names(df[2:end])
        replace!(df[!, n], missing => 0.)
    end
    disallowmissing!(df)
    return df
end
