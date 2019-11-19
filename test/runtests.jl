using Test, ECHOAnalysis, SQLite

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
    datapath = "data"
    db = SQLite.DB(tempname())
    add_taxonomic_profiles(db, datapath, foldermatch="taxprofiles")
    @test "taxa" in SQLite.tables(db).name
    @test SQLite.columns(db, "taxa").name == ["taxon", "abundance", "kind", "sample"]
    @test collect(x.sample for x in SQLite.Query(db, "SELECT DISTINCT sample FROM taxa")) ==  ["C0175_2F_1A","C0192_4F_1A"]

    @test_throws ErrorException add_taxonomic_profiles(db, datapath, foldermatch="taxprofiles")
    # @test_logs (:warn, "removing table taxa") add_taxonomic_profiles(db, datapath, foldermatch="taxprofiles", replace=true)

    add_functional_profiles(db, datapath, kind="genefamilies_relab", foldermatch="funcprofiles")
    @test "genefamilies_relab" in SQLite.tables(db).name
    @test SQLite.columns(db, "genefamilies_relab").name == ["function", "abundance", "stratified", "kind", "sample"]
    @test collect(x.sample for x in SQLite.Query(db, "SELECT DISTINCT sample FROM genefamilies_relab")) ==  ["C0175_2F_1A","C0192_4F_1A"]

    @test_throws ErrorException add_functional_profiles(db, datapath, kind="genefamilies_relab", foldermatch="funcprofiles")
    # @test_logs (:warn, "removing table genefamilies_relab") add_functional_profiles(db, datapath, kind="genefamilies_relab", foldermatch="funcprofiles", replace=true)

    taxa = sqlprofile(db, kind="species")
    @test size(taxa) == (74, 3)
    @test sum(taxa[!, 2]) ≈ 1.
    @test sum(taxa[!, 3]) ≈ 1.

    @test all(!any(ismissing, col) for col in eachcol(taxa))

    func = sqlprofile(db, tablename="genefamilies_relab", kind="genefamilies_relab")
    @test size(func) == (26, 3)
    @test all(!any(ismissing, col) for col in eachcol(func))

    taxa2 = sqlprofile(db, kind="species") do s
        sampleid(s) == "C0175_2F_1A"
    end
    @test size(taxa2)[2] == 2
    @test !any(==(0), taxa2.C0175_2F_1A)
    @test !any(ismissing, taxa2.C0175_2F_1A)
end
