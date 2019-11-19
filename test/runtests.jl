using Test, ECHOAnalysis

sid_strings = ["C0001_1F_1A_S1",
               "C0001_1F_1A",
               "C0002_2F_1A_S2",
               "C0002_1F_1A_S3",
               "C0003_1E_1A_S4",
               "M0004_1E_1A_S5"]
sid_symbols = Symbol.(sid_strings)

@testset "Samples and Timepoints" begin
    @test stoolsample(sid_strings[1]) == StoolSample("C0001_1F_1A", 1, 1, "omnigene")
    @test stoolsample(sid_strings[1]) == stoolsample(sid_strings[2])
    @test stoolsample(sid_symbols[1]) == stoolsample(sid_strings[1])
    @test stoolsample.(sid_strings) == stoolsample.(sid_symbols)

    sids = stoolsample.(sid_strings)

    @test sampleid(sids[2]) == sid_strings[2]
    @test subject(first(sids)) == 1
    @test timepoint(first(sids)) == 1

    @test iskid(first(sids))
    @test ismom(last(sids))
    @test sampletype(first(sids)) == "omnigene"
    @test sampletype(last(sids)) == "ethanol"
    @test_throws ArgumentError stoolsample(replace(sampleid(first(sids)), "_1F_"=>"_1G_"))
    @test_throws ArgumentError stoolsample("C11F1A")

    @test all(x-> x isa AbstractTimepoint, sids)
    @test all(x-> x isa StoolSample, sids)

    sid_letter = resolve_letter_timepoint("42a")
    sids_letter = resolve_letter_timepoint(["42a", "43b", missing])

    @test sid_letter == Timepoint("42a", 42, 1)
    @test resolve_letter_timepoint("42z") == Timepoint("42z", 42, 26)

    @test resolve_letter_timepoint("42A") == sid_letter

    @test sids_letter[1] == sid_letter
    @test ismissing(sids_letter[3])

    @test sortperm(sids) == [1,2,4,3,5,6]
end


@testset "Metadata handling" begin
    sids = stoolsample.(sid_strings)
    @test length(uniquetimepoints(sids)) == 3
    @test length(uniquesubjects(sids)) == 2
    @test length(uniquetimepoints(sids, skipethanol=false)) == 5
    @test length(uniquesubjects(sids, skipethanol=false)) == 4

    @test timepoint(uniquetimepoints(sids)[2]) == 1
    @test timepoint(uniquetimepoints(sids, sortfirst=false)[2]) == 2

    # TODO: Add tests for `breastfeeding` and `formulafeeding`

    @test eltype(numberify(["1", "2", 3])) <: Union{Integer, Missing}
    @test eltype(numberify(["1", "2e3", 3])) <: Union{AbstractFloat, Missing}
    @test eltype(numberify(["1", "2", missing])) <: Union{Integer, Missing}
end

@testset "SQL operations" begin
    # TODO: Write tests for SQL ops
    @test 1==1
end
