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

    ## Profiles

    taxonomic_profiles,
    functional_profiles

    # ## Database Operations
    # sampletable,
    # getlongmetadata,
    # add_taxonomic_profiles,
    # add_functional_profiles,
    # sqlprofile,
    # getallsamples,
    # getmgxmetadata

using DataFrames
using CSV
using Colors
using Microbiome
using BiobakeryUtils
using Dictionaries
using SparseArrays
using ProgressMeter

include("samples.jl")
include("metadata_handling.jl")
include("profiles.jl")
# include("sqlops.jl")

end  # module ECHOAnalysis
