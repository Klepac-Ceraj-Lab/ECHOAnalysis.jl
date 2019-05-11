module ECHOAnalysis

export
    resolve_sampleID,
    add_batch_info,
    merge_tables,
    getmetadatum,
    getmetadata

using DataFrames
using CSV
using PrettyTables

using Microbiome

include("functions.jl")


end  # module ECHOAnalysis
