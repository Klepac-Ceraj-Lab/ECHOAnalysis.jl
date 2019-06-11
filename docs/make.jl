using Documenter

makedocs(
    sitename="ECHO RESONANCE Paper 1",
    pages=[
    "Analysis" => [
        "index.md",
        "01-data-sources.md",
        "02-metadata-merging.md",
        "03-metadata-analysis.md",
        "04-mgx-qc.md",
        "05-mgx-taxonomic-profiles.md",
        "06-omnibus-tests.md"
        ],
    ],
    format = Documenter.HTML(
            prettyurls = get(ENV, "CI", nothing) == "true"
    ),
)
