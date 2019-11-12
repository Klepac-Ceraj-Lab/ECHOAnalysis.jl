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
    numberify

    ## Database Operations

using DataFrames
using MultivariateStats
using CSV
using Colors
using Microbiome

include("samples.jl")
include("metadata_handling.jl")

end  # module ECHOAnalysis
