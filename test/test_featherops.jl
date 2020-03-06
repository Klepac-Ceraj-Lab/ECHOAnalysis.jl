@testset "Feather operations" begin
    look4samples = ["C0175_2F_1A","C0192_4F_1A"]
    datapath = "data"
    testpath = mktempdir(@__DIR__)
    featherall(datapath, testpath, foldermatch="taxprofiles", overwrite=true, skipexisting=true)

    let (indir, outdir) = (readdir(joinpath(datapath * "/taxprofiles")), readdir(testpath))
        @test length(indir) == length(outdir)
        @test map(f-> replace(f, r"\.tsv"=> ".feather"), indir) == outdir
        @test all(f-> sampleid(stoolsample(f)) in look4samples, outdir)
    end


    # @test_throws ErrorException add_taxonomic_profiles(taxdb, datapath, foldermatch="taxprofiles")
    # # @test_logs (:warn, "removing table taxa") add_taxonomic_profiles(db, datapath, foldermatch="taxprofiles", replace=true)
    #
    # funcdb = SQLite.DB()
    # add_functional_profiles(funcdb, datapath, kind="genefamilies_relab", foldermatch="funcprofiles")
    # @test "genefamilies_relab" in SQLite.tables(funcdb).name
    # @test SQLite.columns(funcdb, "genefamilies_relab").name == ["function", "abundance", "stratified", "kind", "sample", "batch"]
    # @test collect(x.sample for x in DBInterface.execute(funcdb, "SELECT DISTINCT sample FROM genefamilies_relab")) ==  look4samples
    #
    # @test_throws ErrorException add_functional_profiles(funcdb, datapath, kind="genefamilies_relab", foldermatch="funcprofiles")
    # # @test_logs (:warn, "removing table genefamilies_relab") add_functional_profiles(db, datapath, kind="genefamilies_relab", foldermatch="funcprofiles", replace=true)
    #
    # taxa = sqlprofile(taxdb, kind="species")
    # @test size(taxa) == (74, 2)
    # @test all(x-> x≈1., sampletotals(taxa))
    #
    # @test all(!ismissing, occurrences(taxa))
    #
    # func = sqlprofile(funcdb, tablename="genefamilies_relab", kind="genefamilies_relab")
    # @test size(func) == (26, 2)
    # @test all(!ismissing, occurrences(func))
    #
    # taxa2 = sqlprofile(taxdb, kind="species") do s
    #     sampleid(s) == look4samples[1]
    # end
    # @test size(taxa2, 2) == 1
    # @test !any(==(0.), occurrences(taxa2))
    # @test !any(ismissing, occurrences(taxa2))
end
