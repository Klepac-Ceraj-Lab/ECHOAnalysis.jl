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
