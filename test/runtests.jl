using Test, ECHOAnalysis, SQLite, Microbiome

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
    taxdb = SQLite.DB()
    funcdb = SQLite.DB()
    add_taxonomic_profiles(taxdb, datapath, foldermatch="taxprofiles")
    @test "taxa" in SQLite.tables(taxdb).name
    @test SQLite.columns(taxdb, "taxa").name == ["taxon", "abundance", "kind", "sample", "batch"]
    @test collect(x.sample for x in SQLite.Query(taxdb, "SELECT DISTINCT sample FROM taxa")) ==  ["C0175_2F_1A","C0192_4F_1A"]

    @test_throws ErrorException add_taxonomic_profiles(taxdb, datapath, foldermatch="taxprofiles")
    # @test_logs (:warn, "removing table taxa") add_taxonomic_profiles(db, datapath, foldermatch="taxprofiles", replace=true)

    add_functional_profiles(funcdb, datapath, kind="genefamilies_relab", foldermatch="funcprofiles")
    @test "genefamilies_relab" in SQLite.tables(funcdb).name
    @test SQLite.columns(funcdb, "genefamilies_relab").name == ["function", "abundance", "stratified", "kind", "sample", "batch"]
    @test collect(x.sample for x in SQLite.Query(funcdb, "SELECT DISTINCT sample FROM genefamilies_relab")) ==  ["C0175_2F_1A","C0192_4F_1A"]

    @test_throws ErrorException add_functional_profiles(funcdb, datapath, kind="genefamilies_relab", foldermatch="funcprofiles")
    # @test_logs (:warn, "removing table genefamilies_relab") add_functional_profiles(db, datapath, kind="genefamilies_relab", foldermatch="funcprofiles", replace=true)

    taxa = sqlprofile(taxdb, kind="species")
    @test size(taxa) == (74, 2)
    @test all(col-> col ≈ 1., sum(occurrences(taxa), dims=1))

    @test all(!any(ismissing, col) for col in eachcol(occurrences(taxa)))

    func = sqlprofile(funcdb, tablename="genefamilies_relab", kind="genefamilies_relab")
    @test size(func) == (26, 2)
    @test all(!any(ismissing, col) for col in eachcol(occurrences(func)))

    taxa2 = sqlprofile(taxdb, kind="species") do s
        sampleid(s) == "C0175_2F_1A"
    end
    @test size(taxa2, 2) == 1
    @test !any(==(0), occurrences(taxa2))
    @test !any(ismissing, occurrences(taxa2))
end
