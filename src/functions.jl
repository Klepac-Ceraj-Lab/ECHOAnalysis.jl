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
                   end
               end
           end
       end
   end
end
