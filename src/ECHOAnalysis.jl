module ECHOAnalysis

export
    ## Samples, Timepoints, and Metdata
    AbstractTimepoint,
    sampleid,
    subject,
    timepoint,
    Timepoint,
    StoolSample,
    stoolsample,
    isstoolsample,
    iskid,
    ismom,
    sampletype,
    replicateid,
    resolve_letter_timepoint,
    uniquetimepoints,
    breastfeeding,
    samplelessthan,
    numberify,
    widemetadata,

    ## Database Operations
    sampletable,
    getlongmetadata,
    add_taxonomic_profiles,
    add_functional_profiles,
    sqlprofile,
    getallsamples,
    getmgxmetadata

using SQLite
using DataFrames
using CSV
using Colors
using Microbiome
using BiobakeryUtils
using Dictionaries
using SparseArrays
using ProgressMeter
using Feather

include("samples.jl")
include("metadata_handling.jl")
include("sqlops.jl")
include("featherops.jl")

end  # module ECHOAnalysis
