using Test, ECHOAnalysis

@testset "Metadata Handling" begin
    sid_strings = ["C0001_1F_1A_S1",
                   "C0001_1F_1A",
                   "C0002_2F_1A_S2",
                   "C0002_1F_1A_S3",
                   "C0003_1E_1A_S4",
                   "M0004_1F_1A_S5"]
    sid_symbols = Symbol.(sid_strings)

    @test resolve_sampleID(sid_strings[1]) == (sample="C0001_1F_1A", subject=1, timepoint=1)
    @test resolve_sampleID(sid_strings[1]) == resolve_sampleID(sid_strings[2])
    @test resolve_sampleID(sid_symbols[1]) == resolve_sampleID(sid_strings[1])
    @test resolve_sampleID(sid_strings) == resolve_sampleID(sid_symbols)

    sids = resolve_sampleID(sid_strings)
    sid = sids[1]

    @test all(x-> typeof(sid) == typeof(x), sids)

    sid_letter = resolve_letter_timepoint("42a")
    sids_letter = resolve_letter_timepoint(["42a", "43b", missing])

    @test sid_letter == (sample="42a", subject=42, timepoint=1)
    @test resolve_letter_timepoint("42z") == (sample="42z", subject=42, timepoint=26)

    @test resolve_letter_timepoint("42A") == sid_letter

    @test sids_letter[1] == sid_letter
    @test ismissing(sids_letter[3])

    @test sortperm(sid_strings, lt=samplelessthan) == [1,2,4,3,5,6]
    @test sortperm(sid_symbols, lt=samplelessthan) == [1,2,4,3,5,6]

    @test length(uniquesamples(sids)) == 4
    @test length(uniquesamples(sids, identifiers=[:subject])) == 3
    @test length(uniquesamples(sids, skipethanol=false)) == 5
    @test length(unique(getfield.(uniquesamples(sids, identifiers=[:subject]), :timepoint))) == 1
    @test length(unique(getfield.(uniquesamples(sids, identifiers=[:subject], takefirst=false), :timepoint))) == 2

    # TODO: Add tests for `breastfeeding` and `formulafeeding`

    @test eltype(numberify(["1", "2", 3])) <: Union{Integer, Missing}
    @test eltype(numberify(["1", "2e3", 3])) <: Union{AbstractFloat, Missing}
    @test eltype(numberify(["1", "2", missing])) <: Union{Integer, Missing}
end
