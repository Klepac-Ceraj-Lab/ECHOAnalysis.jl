module ECHOAnalysis

export
    ## Samples, Timepoints, and Metdata
    AbstractSample,
    Timepoint,
    StoolSample,
    stoolsample,
    uniquesamples,
    breastfeeding,
    samplelessthan,
    numberify,
    widemetadata,
    getsamples

    ## Database Operations

using SQLite
using DataFrames
using CSV
using Colors
using Microbiome

include("samples.jl")
include("metadata_handling.jl")

end  # module ECHOAnalysis
