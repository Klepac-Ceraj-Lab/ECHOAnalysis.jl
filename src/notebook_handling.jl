randrows(df::AbstractDataFrame, nrows=10) = df[rand(nrow(df)) .< nrows / nrow(df), :]
