randrows(df::AbstractDataFrame, nrows=10) = df[rand(nrow(df)) .< nrows / nrow(df), :]

"""
Add notebook information to `data/data.toml`
"""
function updatenotebooks!(tomlpath=tomlpath; notebookspath=joinpath(@__DIR__, "..", "analysis", "notebooks"))
    toml = TOML.parsefile(tomlpath)
    toml["notebooks"] = Dict()
    notebook_files = filter(f-> endswith(f, ".jmd"), readdir(notebookspath))
    for nf in notebook_files
        nbpath = joinpath(notebookspath, nf)

        m = match(r"^(\d+)-", nf)
        isnothing(m) && continue
        n = String(m.captures[1])

        title = ""
        for line in eachline(nbpath)
            !startswith(line, "title:") && continue
            title = replace(replace(line, r"^title:\s?"=>""), "\""=>"")
            break
        end
        nbpath = replace(nbpath, r"^.+\.\.\/"=>"")
        toml["notebooks"][n] = Dict("filepath" => nbpath,
                                    "title"    => title,
                                    "output"   => normpath(joinpath("data", "notebooks", n)),
                                    "figures"  => normpath(joinpath("data", "figures", n)))
    end
    open(tomlpath, "w") do io
        TOML.print(io, toml)
    end
end
