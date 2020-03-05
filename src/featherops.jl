function sparsefeather(datadir::String, filter::Union{Regex,String})
    isdir(datadir) || error("$datadir should be a directory containing taxonmic or functional profiles")
    profiles=String[]
end
