"""
    featherall(source, dest;
            foldermatch="",
            filematch="tsv",
            overwrite=false,
            skipexisting=true)


"""
function featherall(source, dest; foldermatch="", filematch="tsv", overwrite=false, skipexisting=true)
    allfiles = String[]
    for (root, dirs, files) in walkdir(source)
        occursin(foldermatch, root) || continue
        filter!(f-> occursin(filematch, f), files)
        append!(allfiles, joinpath.(root, files))
    end

    @showprogress for f in allfiles
        outfile = joinpath(dest, first(splitext(basename(f))) * ".feather")
        if isfile(outfile)
            skipexisting && continue
            if overwrite
                @warn "overwriting $outfile"
                rm(outfile)
            else
                error("$outfile already exists! Use `overwrite=true` to replace it")
            end
        end 
        Feather.write(outfile, CSV.File(f))
    end
end

function sparsefeather(datadir::String, filter::Union{Regex,String})
    isdir(datadir) || error("$datadir should be a directory containing taxonmic or functional profiles")
    profiles=String[]
end
