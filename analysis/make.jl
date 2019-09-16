using Weave
using Pkg

Pkg.activate(@__DIR__)
inpath = joinpath(@__DIR__, "notebooks")
outpath = joinpath(@__DIR__, "html")

files = filter(f-> endswith(f, ".jmd"), readdir(inpath))
sort!(files)

!isdir(outpath) && mkdir(outpath)
for f in files
    weave(joinpath(inpath, f),
            doctype="github",
            out_path=joinpath(out_path, replace(f, "jmd"=>"html")),
            throw_errors=true
            )
end
