"""
    featherall(source, dest;
            foldermatch="",
            filematch="tsv",
            overwrite=false,
            skipexisting=true)


"""
function featherall(source, dest; foldermatch="", filematch="tsv", overwrite=false, skipexisting=true)
    for (root, dirs, files) in walkdir(source)
        occursin(foldermatch, root) || continue
        filter!(f-> occursin(filematch, f), files)
        for f in files
            infile = joinpath(root, f)
            outfile = joinpath(dest, first(splitext(basename(f))) * ".feather")
            Feather.write(outfile, CSV.File(infile))
        end
    end
end

function sparsefeather(datadir::String, filter::Union{Regex,String})
    isdir(datadir) || error("$datadir should be a directory containing taxonmic or functional profiles")
    profiles=String[]
end
