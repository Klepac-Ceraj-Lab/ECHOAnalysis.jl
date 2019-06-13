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
    getfocusmetadata,
    letter2number,
    parseletterid

using DataFrames
using CSV
using Colors
using ColorBrewer
using PrettyTables

using Microbiome

include("functions.jl")


end  # module ECHOAnalysis
