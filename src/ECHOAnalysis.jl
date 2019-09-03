module ECHOAnalysis

export
    # Metadata handling
    resolve_sampleID,
    resolve_letter_timepoint,
    uniquesamples,
    breastfeeding,
    formulafeeding,
    samplelessthan,
    firstkids,
    numberify,
    getmetadatum,
    getmetadata,
    metacolor,
    getfocusmetadata,

    #File Handling
    add_batch_info,
    merge_tables,

    # Startup Functions
    datatoml,
    color1,
    color2,
    color3,
    color4,
    load_taxonomic_profiles,
    first_kids_tax_profiles,
    load_functional_profiles,
    load_metadata,
    notebookpaths!,

    # Notebook handling
    randrows

using DataFrames
using RecipesBase
using MultivariateStats
using CSV
using Colors
using ColorBrewer
using Pkg.TOML: parsefile
using Microbiome

include_dependency("../data/data.toml")
include("metadata_handling.jl")
include("notebook_handling.jl")
include("file_handling.jl")
include("pcoa_recipe.jl")
include("startup.jl")

end  # module ECHOAnalysis
