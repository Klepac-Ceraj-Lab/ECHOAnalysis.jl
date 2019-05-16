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
        sample = replace(String(m.captures[1]), "-"=>"_")
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
        fname = Symbol(replace(f, suffix=>""))
        rename!(t, names(t)[1]=> :col1, names(t)[2]=>fname)
        df = join(df,t, on=:col1, kind=:outer)
    end
    for n in names(df[2:end])
        df[n] = coalesce.(df[n], 0.)
    end
    disallowmissing!(df)
    return df
end

"""
Get the metadata value for a given subject and timepoint
"""
function getmetadatum(df, metadatum, subject, timepoint=0; default=missing, type=nothing)
    m = filter(df) do row
        !any(ismissing, [row[:metadatum], row[:studyID], row[:timepoint]]) &&
            row[:metadatum] == metadatum &&
            row[:studyID] == subject &&
            row[:timepoint] == timepoint
    end

    nrow(m) == 0 && return default
    nrow(m) > 1 && throw(ErrorException("More than one value found"))

    v = first(m[:value])
    if type === nothing
        return v
    elseif ismissing(v)
        return default
    elseif type <: Number
        return parse(type, v)
    else
        return v
    end
end

"""
Get the metadata values for a given set of subjects and timepoints
"""
function getmetadata(metadf::AbstractDataFrame, subjects::Array{Int,1}, timepoints::Array{Int,1}, metadata::Array{<:AbstractString,1})
    length(subjects) == length(timepoints) || throw(ErrorException("Subjects and timeponts must be the same length"))
    ss = unique(subjects)
    subtps = Dict(s => Set([0, timepoints[findall(isequal(s), subjects)]...]) for s in ss)

    submap = Dict(s => findall(isequal(s), subjects) for s in ss)
    tpmap = Dict(t => findall(isequal(t), timepoints) for t in unique(timepoints))

    metadict = Dict{Int,Dict}()
    notime_metadata = Set()
    timepoint_metadata = Set()

    df = DataFrame(:subject=>subjects, :timepoint=>timepoints)
    for m in Symbol.(metadata)
        df[m] = Array{Any}(missing, length(subjects))
    end

    for row in eachrow(metadf)
        sub = row[:studyID]
        tp = row[:timepoint]
        md = String(row[:metadatum])
        val = row[:value]
        if !any(ismissing, [sub, tp, md]) && md in metadata && sub in subjects && tp in subtps[sub]
            if tp == 0
                rows = submap[sub]
            else
                rows = collect(intersect(submap[sub], tpmap[tp]))
            end
            df[rows, Symbol(md)] = val
        end
    end
    return df
end


"""
Convert metadata into colors
"""
function metacolor(metadata::AbstractVector, colorlevels::AbstractVector; missing_color=colorant"gray")
    um = unique(skipmissing(metadata))
    length(um) <= length(colorlevels) || throw(ErrorException("Not enough colors"))

    color_dict = Dict(um[i] => colorlevels[i] for i in eachindex(um))
    return [ismissing(m) ? missing_color : color_dict[m] for m in metadata]
end
