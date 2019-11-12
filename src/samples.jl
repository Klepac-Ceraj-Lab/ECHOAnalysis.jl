"""
Represents a given timepoint data collection or sample.

Required fields:

- id -- `String`: the unique identifier
- subject -- `Int`: the subject ID
- timepoint -- `Int`: The timepoint.
  For some data types, this will be 0
  indicating that it applies to all timepoints for that subject.
"""
abstract type AbstractTimepoint end

sampleid(s::AbstractTimepoint)  = s.id
subject(s::AbstractTimepoint)   = s.subject
timepoint(s::AbstractTimepoint) = s.timepoint

function Base.show(io::IO, tp::AbstractTimepoint)
    println(io, "Sample id: $(sampleid(tp))")
    println(io, "\tSubject id: $(subject(tp))")
    print(io, "\tTimepoint:  $(timepoint(tp))")
end

struct Timepoint <: AbstractTimepoint
    id::String
    subject::Int
    timepoint::Int
end

"""
Represents a given stool sample

**Fields:**

- id -- `String`: the unique sample identifier
- subject -- `Int`: the subject ID
- timepoint -- `Int`: The timepoint. For some metadata, this will be 0
- type -- `String`: "omnigene" or "ethanol"
"""
struct StoolSample <: AbstractTimepoint
    id::String
    subject::Int
    timepoint::Int
    type::String
end

"""
    stoolsample(sid::Union{AbstractString,Symbol})
    stoolsample(sids::Vector{<:Union{AbstractString,Symbol}})

Parse a sample ID and create a `StoolSample`
or parse a vector of sample IDs and create a vector of `StoolSample`s

Example: For the sample `C0040_3F_1A`:
- `C0040` means ""**C**hild, subject ID **40**"
- `3F` means "timepoint **3**", **F**ecal (omnigene).
  May also be `E` for **E**thanol.
- 1A means "replicate **1**, aliquot **A**"
- Hyphens in sample names are converted to underscores (so "C0001_1F_1A" == "C0001-1F-1A")
- Samples may also end with `_S\\d+` (well #s from sequencing) which are ignored
"""
function stoolsample(sid::String)
    m = match(r"^(([CM])(\d+)[_\-](\d)([FEO])[_\-](\d[A-Z]))(_S\d{1,2})?.?+", sid)

    if m == nothing
        throw(ArgumentError("Input '$sid' doesn't match stool sample"))
    else
        sample = replace(String(m.captures[1]), "-"=>"_")
        subject = parse(Int, m.captures[3])
        timepoint = parse(Int, m.captures[4])
        type = String(m.captures[5])
        if type == "E"
            type = "ethanol"
        elseif type = "F"
            type = "omnigene"
        else
            throw(ArgumentError("unknown type of FecalSample: $type"))
        end

        return StoolSample(id, subject, timepoint, type)
    end
end

stoolsample(sid::Union{AbstractString,Symbol}) = stoolsample(String(sid))
stoolsample(sids::Vector{<:Union{AbstractString,Symbol}}) = stoolsample.(sids)

sampletype(s::StoolSample) = s.type

iskid(s::StoolSample) = occursin(r"^C\d+", s.id)
ismom(s::StoolSample) = occursin(r"^M\d+", s.id)


"""
    resolve_letter_timepoint(sid::AbstractString)
    resolve_letter_timepoint(sid::Missing)

Deals with sampleIDs in which the timepoint is given as a letter.
Eg. `40b` indicates timepoint 2 from subject ID 40. If no letter
is provided, timepoint is assumed to be 1.
"""
function resolve_letter_timepoint(sid::AbstractString)
    m = match(r"^(\d+)([a-zA-Z])?$", sid)
    isnothing(m) && throw(ErrorException("Subject ID has unexpected format: $sid"))
    if isnothing(m.captures[2])
        return (sample = sid,
                subject=parse(Int, sid),
                timepoint=1)
    else
        tp = findfirst(lowercase(m.captures[2]), "abcdefghijklmnopqrstuvwxyz")[1])
        return Timepoint(lowercase(sid), parse(Int, m.captures[1]), tp)
    end
end

resolve_letter_timepoint(sid::Missing) = missing
resolve_letter_timepoint(sids::Vector{<:Union{AbstractString, Missing}}) = resolve_letter_timepoint.(sids)
