"""
    timepointlessthan(x::T, y::T) where T <: AbstractTimepoint
    timepointlessthan(x::T, y::T) where T <: Union{AbstractString, Symbol}

Function to use for sorting `NamedTuple`s from sample IDs.
String versions of sampleIDs will be parsed with `resolve_sampleID`
"""
timepointlessthan(x::T, y::T) where T <: AbstractTimepoint = x.sample < y.sample || (x.sample == y.sample && x.timepoint < y.timepoint)
timepointlessthan(x::T, y::T) where T <: Union{AbstractString, Symbol} = timepointlessthan(stoolsample.((x, y))...)

"""
    uniquesamples(samples::AbstractVector{<:StoolSample}, identifiers::Vector{Symbol}=[:subject,:timepoint];
                        skipethanol=true, samplefilter=x->true)

Identifies unique samples from a vector of sample tuples.
Tuples must have a `sample` field, and additional fields specified in `identifiers`.
Multiple items that have the same values among the `identifiers` will be excluded.

**Example:**

Given the following array of 4 samples:

```@example uniquesamples
s = stoolsample(["C0001_1F_1A", "C0001_1F_2A", "C0001_2F_1A", "C0002_1F_1A"])
```

By default, `identifiers = [:subject, :timepoint]`,
meaning only the second item would be excluded
(since both it and the first item are subject 1, timepoint 1).

A single timepoint for each subject can be acquired by passing `identifiers=[:subject]`.

```@example uniquesamples
uniquesamples(s)
```
```@example uniquesamples
uniquesamples(s, identifiers=[:subject])
```

**Other parameters**

- `skipethanol=true`: exlude samples that match the pattern `_\\dE_`,
    that is ethanol (as opposed to genotek) samples.
- `samplefilter=x->true`: a function to select samples to include.
  By default, all samples are included. Use `samplefilter=iskid`
  to include only child samples for example.
- `takefirst=true`: if `true`, sorts the samples using [`timepointlessthan`]@ref
"""
function uniquesamples(samples::AbstractVector{<:StoolSample};
                        identifiers::Vector{Symbol}=[:subject,:timepoint],
                        skipethanol=true, samplefilter=x->true, takefirst=true)
    seen = NamedTuple[]
    uniquesamples = StoolSample[]
    takefirst && (samples = sort(samples, lt=timepointlessthan))

    map(samples) do s
        !samplefilter(s) && return nothing
        skipethanol && sampletype(s) == "ethanol" && return nothing

        subtp = (; (i => s[i] for i in identifiers)...)
        if !in(subtp, seen)
            push!(seen, subtp)
            push!(uniquesamples, s)
        end
    end
    return uniquesamples
end

function uniquesamples(samples::AbstractVector{<:AbstractString}; kwargs...)
    ss = stoolsample(samples)
    return uniquesamples(ss; kwargs...)
end

"""
    safeoccursin(s::String, thing)
    safeoccursin(r::Regex, thing)

Same as `Base.occursin`, except returns `missing`
if the thing to be searched is `missing`
rather than throwing an error.
"""
safeoccursin(s::String, thing::Any) = occursin(s, thing)
safeoccursin(r::Regex, thing::Any) = occursin(r, thing)
safeoccursin(::String, ::Missing) = missing
safeoccursin(::Regex, ::Missing) = missing

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
function breastfeeding(row::DataFrameRow)
    bf = any([
        safeoccursin(r"[Bb]reast", row[:milkFeedingMethods]),
        safeoccursin(r"[Yy]es", row[:exclusivelyNursed]),
        all(map(x->safeoccursin(r"[Nn]o", x), row[[:exclusiveFormulaFed, :exclusivelyNursed]])),
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
function formulafeeding(row::DataFrameRow)
    ff = any([
        safeoccursin(r"[Ff]ormula", row[:milkFeedingMethods]),
        safeoccursin(r"[Yy]es", row[:exclusiveFormulaFed]),
        all(map(x->safeoccursin(r"[Nn]o",x), row[[:exclusiveFormulaFed, :exclusivelyNursed]])),
        !ismissing(row[:amountFormulaPerFeed]),
        !ismissing(row[:formulaTypicalType]),
        row[:amountFormulaPerFeed] > 0
        ])
    return !ismissing(ff) && ff
end


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
        safeoccursin(r"[\.e]", x) ? parse(Float64, x) : parse(Int, x)
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
        !any(ismissing, [row[:metadatum], row[:subject], row[:timepoint]]) &&
            row[:metadatum] == metadatum &&
            row[:subject] == subject &&
            row[:timepoint] == timepoint
    end

    nrow(m) == 0 && return default
    nrow(m) > 1 && throw(ErrorException("More than one value found"))

    v = first(m[!,:value])
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
        df[!,m] = Array{Any}(missing, length(subjects))
    end

    for row in eachrow(metadf)
        sub = row[:subject]
        tp = row[:timepoint]
        md = String(row[:metadatum])
        val = row[:value]
        if !any(ismissing, [sub, tp, md]) && md in metadata && sub in subjects && tp in subtps[sub]
            if tp == 0
                rows = submap[sub]
            else
                rows = collect(intersect(submap[sub], tpmap[tp]))
            end
            df[rows, Symbol(md)] .= val
        end
    end
    return df
end

const metadata_focus_headers = String[
    # fecal sample info
    "Mgx_batch",
    "DOM",
    "RAInitials_Extract",
    # basic data
    "childGender",
    "milkFeedingMethods",
    "APOE",
    "birthType",
    "exclusivelyNursed",
    "exclusiveFormulaFed",
    "formulaTypicalType",
    "mother_HHS",
    "motherHHS_Occu",
    "motherHHS_Edu",
    "assessmentDate",
    # calculated
    "ageLabel",
    "cogAssessment",
    "cogScore",
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
    "subcortical",
    "neocortical",
    "cerebellar",
    "limbic",
    "hires_total",
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
    ### WISC (6-16 years)
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
customprocess(col, ::MDColumn{:cerebellar})                         = numberify(col)
customprocess(col, ::MDColumn{:subcortical})                        = numberify(col)
customprocess(col, ::MDColumn{:neocortical})                        = numberify(col)
customprocess(col, ::MDColumn{:limbic})                             = numberify(col)
customprocess(col, ::MDColumn{:hires_total})                        = numberify(col)
customprocess(col, ::MDColumn{:white_matter_volume})                = numberify(col)
customprocess(col, ::MDColumn{:grey_matter_volume})                 = numberify(col)
customprocess(col, ::MDColumn{:csf_volume})                         = numberify(col)


## Cog scores
customprocess(col, ::MDColumn{:cogScore})                           = numberify(col)
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
function customprocess(col, ::MDColumn{:mother_HHS})
    col = numberify(col)
    # some missing entries are encoded as 9999
    for i in eachindex(col)
        if !ismissing(col[i]) && col[i] >= 9999
            col[i] = missing
        end
    end
    return col
end

function customprocess(col, ::Union{MDColumn{:motherHHS_Occu},MDColumn{:motherHHS_Edu}})
    return customprocess(col, MDColumn(:mother_HHS))
end
