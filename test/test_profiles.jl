@testset "Profiles" begin
    look4samples = ["C0175_2F_1A","C0192_4F_1A"]
    datapath = "data"

    taxa = taxonomic_profiles(joinpath(datapath, "taxprofiles"))
    @test size(taxa) == (74, 2)
    @test all(x-> x≈1., sampletotals(taxa))

    @test all(!ismissing, occurrences(taxa))

    func = functional_profiles(joinpath(datapath, "funcprofiles"), kind="genefamilies_relab")
    @test size(func) == (26, 2)
    @test all(!ismissing, occurrences(func))

    taxa2 = sqlprofile(taxdb, kind="species") do s
        sampleid(s) == look4samples[1]
    end

    @test size(taxa2, 2) == 1
    @test !any(==(0.), occurrences(taxa2))
    @test !any(ismissing, occurrences(taxa2))
end
