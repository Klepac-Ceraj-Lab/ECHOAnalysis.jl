using Test, ECHOAnalysis, SQLite, Microbiome

sid_strings = ["C0001_1F_2A",
               "C0001_1F_1A_S1",
               "C0001_1F_1A",
               "C0002_2F_1A_S2",
               "C0002_1F_1A_S3",
               "C0003_1E_1A_S4",
               "M0004_1E_1A_S5"]
sid_symbols = Symbol.(sid_strings)

include("test_samples.jl")
include("test_metadata_handling.jl")
include("test_sqlops.jl")
include("test_featherops.jl")
