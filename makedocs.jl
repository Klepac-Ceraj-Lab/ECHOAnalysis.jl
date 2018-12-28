using Documenter

makedocs(
    sitename="ECHO RESONANCE Paper 1",
    pages=[
    "Analysis" => [
        "About this section" => "notebooks/index.md",
        "Data description" => "notebooks/01-data-sources.md",
        "Metadata Summary" => "notebooks/02-metadata-summary.md",
        ],
    ],
    Documenter.HTML(
        prettyurls = false
    )
)
