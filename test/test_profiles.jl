@testset "Profiles" begin
    look4samples = ["C0175_2F_1A","C0192_4F_1A"]
    datapath = "data"

    longtax, taxfeatures, taxsamples = taxonomic_profiles(joinpath(datapath, "taxprofiles"))
    @test length(taxfeatures) == 74
    @test length(taxsamples) == 2
    @test Set(longtax.taxon) == taxfeatures
    @test Set(longtax.sample) == taxsamples

    taxa = widen2comm(longtax, taxfeatures, taxsamples)
    @test size(taxa) == (74, 2)
    @test_broken all(x-> x≈1., sampletotals(taxa))
    @test all(!ismissing, occurrences(taxa))

    longfunc, funcfeatures, funcsamples = functional_profiles(joinpath(datapath, "funcprofiles"), kind="genefamilies_relab")
    @test length(funcfeatures) == 26
    @test length(funcsamples) == 2
    @test Set(longfunc.func) == funcfeatures
    @test Set(longfunc.sample) == funcsamples

    func = widen2comm(longfunc, funcfeatures, funcsamples, featurecol=:func)
    @test size(func) == (26, 2)
    @test all(!ismissing, occurrences(func))
end
