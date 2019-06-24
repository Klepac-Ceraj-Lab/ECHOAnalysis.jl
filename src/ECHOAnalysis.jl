module ECHOAnalysis

export
    # Metadata handling
    resolve_sampleID,
    resolve_letter_timepoint,
    breastfeeding,
    formulafeeding,
    firstkids,
    numberify,

    getmetadatum,
    getmetadata,
    metacolor,
    convert2num,
    getfocusmetadata,
    letter2number,
    parseletterid,

    #File Handling
    add_batch_info,
    merge_tables,


using DataFrames
using CSV
using Colors
using ColorBrewer
using PrettyTables

using Microbiome

include("metadata_handling.jl")
include("file_handling.jl")

end  # module ECHOAnalysis
