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

"""
    sampleid(s::AbstractTimepoint)

Get the `id` field from an `AbstractTimepoint`.
"""
sampleid(s::AbstractTimepoint)  = s.id

"""
    subject(s::AbstractTimepoint)

Get the `subject` field from an `AbstractTimepoint`.
"""
subject(s::AbstractTimepoint)   = s.subject

"""
    timepoint(s::AbstractTimepoint)

Get the `timepoint` field from an `AbstractTimepoint`.
"""
timepoint(s::AbstractTimepoint) = s.timepoint

function Base.isless(x::AbstractTimepoint, y::AbstractTimepoint)
    subject(x) < subject(y) || (subject(x) == subject(y) && timepoint(x) < timepoint(y))
end

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
- replicate -- `String`: replicate / aliquot ID (eg `"2B"`)
- type -- `String`: "omnigene" or "ethanol"
"""
struct StoolSample <: AbstractTimepoint
    id::String
    subject::Int
    timepoint::Int
    replicate::String
    type::String
end

"""
    isstoolsample(sid::Union{AbstractString, AbstractTimepoint, StoolSample})

Returns `true` if sid is a `StoolSample`,
or if an `AbstractTimepoint`'s `id` field or `AbstractString`
conforms to the expected form of an ECHO RESONANCE stool sample.
"""
isstoolsample(::StoolSample) = true
isstoolsample(sid::AbstractTimepoint) = occursin(r"^(([CM])(\d+)[_\-](\d+)([FEO])[_\-](\d[A-Z]))(_S\d{1,2})?.?+", sampleid(sid))
isstoolsample(sid::AbstractString) = occursin(r"^(([CM])(\d+)[_\-](\d+)([FEO])[_\-](\d[A-Z]))(_S\d{1,2})?.?+", sid)

"""
    stoolsample(sid::Union{AbstractString,Symbol})

Parse a sample ID and create a `StoolSample`.

Example: For the sample `C0040_3F_1A`:
- `C0040` means ""**C**hild, subject ID **40**"
- `3F` means "timepoint **3**", **F**ecal (omnigene).
  May also be `E` for **E**thanol.
- 1A means "replicate **1**, aliquot **A**"
- Hyphens in sample names are converted to underscores (so "C0001_1F_1A" == "C0001-1F-1A")
- Samples may also end with `_S\\d+` (well #s from sequencing) which are ignored
"""
function stoolsample(sid::String)
    isstoolsample(sid) || throw(ArgumentError("Input '$sid' doesn't match stool sample"))
    m = match(r"^(([CM])(\d+)[_\-](\d+)([FEO])[_\-](\d[A-Z]))(_S\d{1,2})?.?+", sid)

    sample = replace(String(m.captures[1]), "-"=>"_")
    subject = parse(Int, m.captures[3])
    timepoint = parse(Int, m.captures[4])
    type = String(m.captures[5])
    replicate = String(m.captures[6])
    if type == "E"
        type = "ethanol"
    elseif type == "F"
        type = "omnigene"
    else
        throw(ArgumentError("unknown type of FecalSample: $type"))
    end

    return StoolSample(sample, subject, timepoint, replicate, type)
end

stoolsample(sid::Union{AbstractString,Symbol}) = stoolsample(String(sid))

"""
    sampletype(s::Union{AbstractString, StoolSample})

Get the `type` field from a `StoolSample`
(or a `String` that can be converted to `StoolSample`).
"""
sampletype(s::StoolSample) = s.type
sampletype(s::AbstractString) = sampletype(stoolsample(s))

"""
    replicateid(s::Union{AbstractString, StoolSample})

Get the `replicate` field from a `StoolSample`
(or a `String` that can be converted to `StoolSample`).
"""
replicateid(s::StoolSample) = s.replicate
replicateid(s::AbstractString) = replicateid(stoolsample(s))

"""
    iskid(s::Union{AbstractString, StoolSample})

Return `true` if `StoolSample`
(or a `String` that can be converted to `StoolSample`)
refers to a child sample, `false` otherwise.
"""
iskid(s::StoolSample) = occursin(r"^C\d+", s.id)
iskid(s::AbstractString) = iskid(stoolsample(s))

"""
    ismom(s::Union{AbstractString, StoolSample})

Return `true` if `StoolSample`
(or a `String` that can be converted to `StoolSample`)
refers to a maternal sample, `false` otherwise.
"""
ismom(s::StoolSample) = occursin(r"^M\d+", s.id)
ismom(s::AbstractString) = ismom(stoolsample(s))

function Base.isless(x::StoolSample, y::StoolSample)
    if subject(x) < subject(y)
        return true
    elseif subject(x) == subject(y)
        if timepoint(x) < timepoint(y)
            return true
        elseif timepoint(x) == timepoint(y)
            return replicateid(x) < replicateid(y)
        else
            return false
        end
    else
        return false
    end
end

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
        return Timepoint(string(sid),parse(Int, sid), 1)
    else
        tp = findfirst(lowercase(m.captures[2]), "abcdefghijklmnopqrstuvwxyz")[1]
        return Timepoint(lowercase(sid), parse(Int, m.captures[1]), tp)
    end
end

resolve_letter_timepoint(sid::Missing) = missing
