"""
    uniquetimepoints(samples::AbstractVector{<:StoolSample};
                        skipethanol=true,
                        samplefilter=x->true,
                        sortfirst=true,
                        takefirst=true)

Identifies unique timepoints (that is, `subject=>timepoint` pairs)
from a vector of [`AbstractTimepoint`]@ref.

**Example:**

Given the following array of 4 samples:

```@example uniquetimepoints
s = stoolsample(["C0001_1F_1A", "C0001_1F_2A", "C0001_2F_1A", "C0002_1F_1A"])
```

The first three samples are from the same `subject`,
the first two of which are from the same timepoint (they're replicates).

By default, this function excludes duplicate subject-timepoints,
but keeps multiple timepoints for each individual subject.
If a single timepoint from each subject is desired,
use `takefirst=true`.
In this later case, leave `sortfirst=true` (the default)
to be sure to take the earliest timepoint for each subject.

**Other parameters**

- `skipethanol=true`: exlude samples that match the pattern `_\\dE_`,
    that is ethanol (as opposed to genotek) samples.
- `samplefilter=x->true`: a function to select samples to include.
  By default, all samples are included. Use `samplefilter=iskid`
  to include only child samples for example.
- `sortfirst=true`: if `true`, sorts the vector at the beginning
- `takefirst=false`: if `true`, only takes 1 timepoint for each subject.
    Use with `sortfirst` to get only the first timepoint.
"""
function uniquetimepoints(samples::AbstractVector{<:StoolSample};
                           skipethanol=true, samplefilter=x->true,
                           sortfirst=true, takefirst=false)
    seen = Tuple[]
    uniquesamples = StoolSample[]
    sortfirst && (samples = sort(samples))

    map(samples) do s
        !samplefilter(s) && return nothing
        skipethanol && sampletype(s) == "ethanol" && return nothing

        takefirst ? subtp = (subject(s),) : subtp = (subject(s), timepoint(s))
        if !in(subtp, seen)
            push!(seen, subtp)
            push!(uniquesamples, s)
        end
    end
    return uniquesamples
end

function uniquetimepoints(samples::AbstractVector{<:AbstractString}; kwargs...)
    ss = stoolsample.(samples)
    return uniquetimepoints(ss; kwargs...)
end

function airtable_request(method, cred, path; query_kwargs...)
    query = ["api_key"=>cred]
    for (key, value) in query_kwargs
        isempty(value) && continue
        push!(query, string(key) => string(value))
    end
    uri = HTTP.URI(host="api.airtable.com", scheme="https", path=path, query=query)
    resp = HTTP.request(method, uri)
    return JSON3.read(String(resp.body))
end

"""
    airtable_metadata(key=ENV["AIRTABLE_KEY"])

Get fecal sample metadata table from airtable.

The API `key` comes from https://airtable.com/account.
"""
function airtable_metadata(key=ENV["AIRTABLE_KEY"])
    records = []
    req = airtable_request("GET", key, "/v0/appyRaPsZ5RsY4A1h/Master"; view="Everything", filterByFormula="NOT({Mgx_batch}='')")
    append!(records, req.records)
    while haskey(req, :offset) && length(records) < 2200
        @info "Making another request"
        req = airtable_request("GET", key, "/v0/appyRaPsZ5RsY4A1h/Master"; view="Everything", filterByFormula="NOT({Mgx_batch}='')", offset=req.offset)
        append!(records, req.records)
        sleep(0.250)
    end
    vcat(map(records) do record
        f = record.fields
        DataFrame(sample    = f.SampleID,
                  subject   = parse(Int, f.SubjectID),
                  timepoint = parse(Int, f.TimePoint),
                  batch     = :Mgx_batch in keys(f) ? parse(Int, match(r"Batch (\d+)", f.Mgx_batch).captures[1]) : missing
                  )
    end...)
end
