using Documenter

makedocs(
    sitename="ECHO RESONANCE Paper 1",
    pages=[
    "Analysis" => [
        "About this section" => "index.md",
        "Data description" => "01-data-sources.md",
        "Metadata Merging" => "02-metadata-merging.md",
        # "Metadata Analysis" => "03-metadata-analysis.md",
        ],
    ],
    format = Documenter.HTML(
            prettyurls = get(ENV, "CI", nothing) == "true"
    )
)
