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
        Feather.write(outfile, CSV.File(f))
    end
    end
end

function sparsefeather(datadir::String, filter::Union{Regex,String})
    isdir(datadir) || error("$datadir should be a directory containing taxonmic or functional profiles")
    profiles=String[]
end
