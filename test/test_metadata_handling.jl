@testset "Metadata handling" begin
    sids = stoolsample.(sid_strings)
    @test length(uniquetimepoints(sids)) == 3
    @test length(uniquetimepoints(sids, takefirst=true)) == 2
    @test length(uniquetimepoints(sids, skipethanol=false)) == 5
    @test length(uniquetimepoints(sids, takefirst=true, skipethanol=false)) == 4

    @test timepoint(uniquetimepoints(sids)[2]) == 1
    @test timepoint(uniquetimepoints(sids, sortfirst=false)[2]) == 2

    # TODO: Add tests for `breastfeeding` and `formulafeeding`

    @test eltype(numberify(["1", "2", 3])) <: Union{Integer, Missing}
    @test eltype(numberify(["1", "2e3", 3])) <: Union{AbstractFloat, Missing}
    @test eltype(numberify(["1", "2", missing])) <: Union{Integer, Missing}
end
