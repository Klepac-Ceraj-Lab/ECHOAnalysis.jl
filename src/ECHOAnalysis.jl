module ECHOAnalysis

export
    resolve_sampleID,
    add_batch_info,
    merge_tables,
    getmetadatum,
    getmetadata,
    metacolor,
    convert2num,
    breastfeeding,
    formulafeeding,
    firstkids,
    numberify,
    getfocusmetadata

using DataFrames
using CSV
using Colors
using ColorBrewer
using PrettyTables

using Microbiome

include("functions.jl")


end  # module ECHOAnalysis
