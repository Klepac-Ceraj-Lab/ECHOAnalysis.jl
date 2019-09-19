using Weave
using Pkg

Pkg.activate(@__DIR__)
inpath = joinpath(@__DIR__, "notebooks")
outpath = joinpath(@__DIR__, "html")

files = filter(f-> endswith(f, ".jmd"), readdir(inpath))
sort!(files)

!isdir(outpath) && mkdir(outpath)
for f in sort(files)[1:9]
    @info "Weaving $f"
    weave(joinpath(inpath, f),
            out_path=joinpath(outpath, replace(f, "jmd"=>"html")),
            throw_errors=true
            )
end
