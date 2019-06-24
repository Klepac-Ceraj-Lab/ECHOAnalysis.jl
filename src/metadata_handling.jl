"""
    resolve_sampleID(sid::Union{AbstractString,Symbol})

Parse a sample name and return its components as a named tuple:

Example: For the sample `C0040_3F_1A`:
- `C0040` means ""**C**hild, subject ID **40**"
- `3F` means "timepoint **3**", **F**ecal (genotech).
  May also be `E` for **E**thanol.
- 1A means "replicate **1**, aliquot **A**"

This function returns the `NamedTuple` `(sample="C0040_3F_1A", subject=40, timepoint=3)`
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
    resolve_letter_timepoint(sid::AbstractString)
    resolve_letter_timepoint(sid::Missing)

Deals with sampleIDs in which the timepoint is given as a letter.
Eg. `40b` indicates timepoint 2 from subject ID 40. If no letter
is provided, timepoint is assumed to be 1.
"""
function resolve_letter_timepoint(sid::AbstractString)
    m = match(r"(\d+)([a-zA-Z])?", sid)
    isnothing(m) && throw(ErrorException("Subject ID has unexpected format: $sid"))
    if isnothing(m.captures[2])
        return (sample = s,
                subject=parse(Int, sid)
                timepoint=1)
    else
        return (sample  = s,
                subject = parse(Int, m.captures[1]),
                timepoint = findfirst(
                                lowercase(m.captures[2]),
                                "abcdefghijklmnopqrstuvwxyz")[1])
    end
end
resolve_letter_timepoint(sid::Missing) = missing

import Base.occursin

# WARNING: This is type piracy.
# But it beats having to check if the thing is missing each time
occursin(::String, ::Missing) = missing
occursin(::Regex, ::Missing) = missing

"""
    breastfeeding(row::DataFrameRow)

Checks a wide-form metadata row for breastfeeding information.
Returns `true` if any of the following are `true`:
- `:milkFeedingMethods` contains "breast"
- `:exclusivelyNursed` is "Yes" (or "yes")
- Both `:exclusivelyNursed` and `:exclusiveFormulaFed` are "No" (or "no")
- Any of the following have values > 0:
    - `:typicalNumberOfFeedsFromBreast`
    - `:typicalNumberOfEpressedMilkFeeds`
    - `:lengthExclusivelyNursedMonths`
    - `:noLongerFeedBreastmilkAge`

Otherwise returns `false`.
"""
function breastfeeding(row)
    bf = any([
        occursin(r"[Bb]reast", row[:milkFeedingMethods]),
        occursin(r"[Yy]es", row[:exclusivelyNursed]),
        all(occursin.(r"[Nn]o", row[[:exclusiveFormulaFed, :exclusivelyNursed]])),
        row[:typicalNumberOfFeedsFromBreast] > 0,
        row[:typicalNumberOfEpressedMilkFeeds] > 0,
        row[:lengthExclusivelyNursedMonths] > 0,
        row[:noLongerFeedBreastmilkAge] > 0
        ])
    return !ismissing(bf) && bf
end

"""
    formulafeeding(row::DataFrameRow)

Checks a wide-form metadata row for formula feeding information.
Returns `true` if any of the following are `true`:
- `:milkFeedingMethods` contains "Formula" (or "formula")
- `:exclusivelyFormulaFed` is "Yes" (or "yes")
- Both `:exclusivelyNursed` and `:exclusiveFormulaFed` are "No" (or "no")
- `:amountFormulaPerFeed` is > 0:
- Either `:amountFormulaPerFeed` or `:formulaTypicalType` have (non-missing) values

Otherwise returns `false`.
"""
function formulafeeding(row)
    ff = any([
        occursin(r"[Ff]ormula", row[:milkFeedingMethods]),
        occursin(r"[Yy]es", row[:exclusiveFormulaFed]),
        all(occursin.(r"[Nn]o", row[[:exclusiveFormulaFed, :exclusivelyNursed]])),
        !ismissing(row[:amountFormulaPerFeed]),
        !ismissing(row[:formulaTypicalType]),
        row[:amountFormulaPerFeed] > 0
        ])
    return !ismissing(ff) && ff
end

"""
    firstkids(samples::Vector{<:NamedTuple})
    firstkids(samples::Vector{<:AbstractString})

From a list of sample ids, identify the earliest sample for each child.
`Samples` may be sample IDs that can be parsed by [`resolve_sampleID`](@ref),
or `NamedTuple`s containing `:subject` and `:timepoint` fields.

Returns a vector of indicies that can be used to slice the original vector.
"""
function firstkids(samples::AbstractVector{<:NamedTuple})
    subjects = Dict()

    for i in eachindex(samples)
        s = samples[i]
        startswith(s[:sample], "C") || continue
        if s[:subject] in keys(subjects)
            if s[:timepoint] < subjects[s[:subject]][:timepoint]
                continue
            end
        end

        subjects[s[:subject]] = (timepoint=s[:timepoint], index=i)
    end
    [subjects[k][:index] for k in keys(subjects)]
end

firstkids(samples::AbstractVector{<:AbstractString}) = firstkids(resolve_sampleID.(samples))


"""
    numberify(x)
    numberify(v::AbstractArray)

Try to convert something into a number or array of numbers using the heuristics:
- if x is already a number or is `missing`, return x
- if v is already an array of numbers, return v
- if x is a `String`:
    - if x contains a `.` or `e`, parse as a `Float64`
    - Otherwise parse as an Int
- if numberified vector contains any `Float`s, return a vector of `missing` and `Float64`
- if x is anything other than a String or a number, return `missing`
"""
function numberify(x)
    ismissing(x) && return missing
    x isa Real && return x
    if x isa AbstractString
        occursin(r"[\.e]", x) ? parse(Float64, x) : parse(Int, x)
    else
        @warn "Something weird" x typeof(x)
        return missing
    end
end

function numberify(v::AbstractArray)
    eltype(v) <: Union{Missing, <:Real} && return v
    v = numberify.(v)
    if any(i -> i isa Float64, v)
        return Union{Missing, Float64}[y for y in v]
    else
        return Union{Missing, Int}[y for y in v]
    end
end

"""
    metacolor(metadata::AbstractVector, colorlevels::AbstractVector; missing_color=colorant"gray")

Convert metadata into colors
"""
function metacolor(metadata::AbstractVector, colorlevels::AbstractVector; missing_color=colorant"gray")
    um = unique(skipmissing(metadata))
    length(um) <= length(colorlevels) || throw(ErrorException("Not enough colors"))

    color_dict = Dict(um[i] => colorlevels[i] for i in eachindex(um))
    return [ismissing(m) ? missing_color : color_dict[m] for m in metadata]
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

const metadata_focus_headers = String[
    "childGender",
    "milkFeedingMethods",
    "APOE",
    "birthType",
    "exclusivelyNursed",
    "exclusiveFormulaFed",
    "formulaTypicalType",
    "motherSES",
    "assessmentDate",
    # numeric
    "correctedAgeDays",
    "amountFormulaPerFeed",
    "lengthExclusivelyNursedMonths",
    "typicalNumberOfEpressedMilkFeeds",
    "typicalNumberOfFeedsFromBreast",
    "noLongerFeedBreastmilkAge",
    "ageStartSolidFoodMonths",
    "childHeight",
    "childWeight",
    ## From brain data
    "LThickness",
    "RThickness",
    "LSurfArea",
    "RSurfArea",
    "ICV",
    "subcortical_volume",
    "white_matter_volume",
    "grey_matter_volume",
    "csf_volume",
    ## Cognitive assessments
    ### Mullen (0-3 years 11 months)
    "mullen_VerbalComposite",
    "mullen_NonVerbalComposite",
    "mullen_EarlyLearningComposite", # average of other two
    ### Bayley's (1-25 months)
    "adaptiveBehaviorComposite",
    "languageComposite",
    "socialEmotionalComposite",
    "motorComposite",
    ### WISC- (6-16 years)
    "FRI_CompositeScore",
    "PSI_Composite",
    "VCI_CompositeScore",
    "VSI_CompositeScore",
    "WMI_Composite",
    "FSIQ_Composite", # average
    ### WPPSI-IV (4-5 years 11 months)
    "fluidReasoningComposite",
    "verbalComprehensionComposite",
    "visualSpatialComposite",
    "workingMemoryComposite",
    "fullScaleComposite", # average
    ]

# Use value types for special cases. See
#   https://docs.julialang.org/en/v1/manual/types/#%22Value-types%22-1
#   https://discourse.julialang.org/t/special-case-handling-alternatives-to-if-else/21608/4
struct MDColumn{T}
end

MDColumn(s::String) = MDColumn{Symbol(s)}()
MDColumn(s::Symbol) = MDColumn{s}()

customprocess(col, ::MDColumn) = col

# Make numeric
customprocess(col, ::MDColumn{:correctedAgeDays})                   = numberify(col)
customprocess(col, ::MDColumn{:amountFormulaPerFeed})               = numberify(col)
customprocess(col, ::MDColumn{:lengthExclusivelyNursedMonths})      = numberify(col)
customprocess(col, ::MDColumn{:typicalNumberOfEpressedMilkFeeds})   = numberify(col)
customprocess(col, ::MDColumn{:typicalNumberOfFeedsFromBreast})     = numberify(col)
customprocess(col, ::MDColumn{:noLongerFeedBreastmilkAge})          = numberify(col)
customprocess(col, ::MDColumn{:ageStartSolidFoodMonths})            = numberify(col)
customprocess(col, ::MDColumn{:childHeight})                        = numberify(col)

## brain volumes
customprocess(col, ::MDColumn{:LThickness})                         = numberify(col)
customprocess(col, ::MDColumn{:RThickness})                         = numberify(col)
customprocess(col, ::MDColumn{:LSurfArea})                          = numberify(col)
customprocess(col, ::MDColumn{:RSurfArea})                          = numberify(col)
customprocess(col, ::MDColumn{:ICV})                                = numberify(col)
customprocess(col, ::MDColumn{:subcortical_volume})                 = numberify(col)
customprocess(col, ::MDColumn{:white_matter_volume})                = numberify(col)
customprocess(col, ::MDColumn{:grey_matter_volume})                 = numberify(col)
customprocess(col, ::MDColumn{:csf_volume})                         = numberify(col)




### Mullen
customprocess(col, ::MDColumn{:mullen_VerbalComposite})             = numberify(col)
customprocess(col, ::MDColumn{:mullen_NonVerbalComposite})          = numberify(col)
customprocess(col, ::MDColumn{:mullen_EarlyLearningComposite})      = numberify(col)
### WSSPI-IV
customprocess(col, ::MDColumn{:fluidReasoningComposite})            = numberify(col)
customprocess(col, ::MDColumn{:fullScaleComposite})                 = numberify(col)
customprocess(col, ::MDColumn{:verbalComprehensionComposite})       = numberify(col)
customprocess(col, ::MDColumn{:visualSpatialComposite})             = numberify(col)
customprocess(col, ::MDColumn{:workingMemoryComposite})             = numberify(col)
### WISC-V
customprocess(col, ::MDColumn{:FRI_CompositeScore})                 = numberify(col)
customprocess(col, ::MDColumn{:FSIQ_Composite})                     = numberify(col)
customprocess(col, ::MDColumn{:PSI_Composite})                      = numberify(col)
customprocess(col, ::MDColumn{:VCI_CompositeScore})                 = numberify(col)
customprocess(col, ::MDColumn{:VSI_CompositeScore})                 = numberify(col)
customprocess(col, ::MDColumn{:WMI_Composite})                      = numberify(col)
### Bayley's
customprocess(col, ::MDColumn{:adaptiveBehaviorComposite})          = numberify(col)
customprocess(col, ::MDColumn{:languageComposite})                  = numberify(col)
customprocess(col, ::MDColumn{:socialEmotionalComposite})           = numberify(col)
customprocess(col, ::MDColumn{:motorComposite})                     = numberify(col)

# other custom processing
function customprocess(col, ::MDColumn{:motherSES})
    col = numberify(col)
    # some missing entries are encoded as 9999
    for i in eachindex(col)
        if !ismissing(col[i]) && col[i] == 9999
            col[i] = missing
        end
    end
    return col
end

"""
    getfocusmetadata(df, samples; focus=metadata_focus_headers)

Get a wide-form metadata table with a subset of headers for a subset of subject/sample IDs
"""
function getfocusmetadata(df::AbstractDataFrame, samples::Vector{<:NamedTuple}; focus=metadata_focus_headers)
    subjects = [s.subject for s in samples]
    timepoints = [s.timepoint for s in samples]
    @show first(subjects)
    @show first(timepoints)
    df = getmetadata(df, subjects, timepoints, metadata_focus_headers)

    for n in names(df)
        df[n] = customprocess(df[n], MDColumn(n))
    end

    return df
end


"""
    getfocusmetadata(longfilepath, samples; focus=metadata_focus_headers)

Get a wide-form metadata table with a subset of headers for a subset of subject/sample IDs
"""
function getfocusmetadata(longfilepath::AbstractString, samples::Vector{<:NamedTuple}; focus=metadata_focus_headers)
    md = CSV.read(longfilepath)
    getfocusmetadata(md, samples; focus=focus)
end

#
function getfocusmetadata(longfilepath, samples::Vector{<:AbstractString}; focus=metadata_focus_headers)
    getfocusmetadata(longfilepath, resolve_sampleID.(samples), focus=focus)
end
