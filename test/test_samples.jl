@testset "Samples and Timepoints" begin
    @test stoolsample(sid_strings[1]) == StoolSample("C0001_1F_2A", 1, 1, "2A", "omnigene")
    @test stoolsample(sid_strings[2]) == stoolsample(sid_strings[3])
    @test stoolsample(sid_symbols[1]) == stoolsample(sid_strings[1])
    @test stoolsample.(sid_strings) == stoolsample.(sid_symbols)

    sids = stoolsample.(sid_strings)

    @test sampleid(sids[3]) == sid_strings[3]
    @test subject(first(sids)) == 1
    @test timepoint(first(sids)) == 1
    @test replicateid(first(sids)) == "2A"

    @test iskid(first(sids))
    @test ismom(last(sids))
    @test sampletype(first(sids)) == "omnigene"
    @test sampletype(last(sids)) == "ethanol"
    @test_throws ArgumentError stoolsample(replace(sampleid(first(sids)), "_1F_"=>"_1G_"))
    @test_throws ArgumentError stoolsample("C11F1A")

    @test all(x-> x isa AbstractTimepoint, sids)
    @test all(x-> x isa StoolSample, sids)

    sid_letter = resolve_letter_timepoint("42a")
    sids_letter = resolve_letter_timepoint.(["42a", "43b", missing])

    @test sid_letter == Timepoint("42a", 42, 1)
    @test resolve_letter_timepoint("42z") == Timepoint("42z", 42, 26)

    @test resolve_letter_timepoint("42A") == sid_letter

    @test sids_letter[1] == sid_letter
    @test ismissing(sids_letter[3])

    @test sortperm(sids) == [2,3,1,5,4,6,7]

    @test stoolsample("C0001_1F_1A") < stoolsample("C0001_1F_2A")
end
