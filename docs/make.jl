using Documenter, ECHOAnalysis

makedocs(
    sitename="ECHO RESONANCE Paper 1",
    pages=[
    "Getting Started" => "index.md",
    "Dealing with Metadata" => "metadata_handling.md",
        ],
    format = Documenter.HTML(
            prettyurls = get(ENV, "CI", nothing) == "true"
    ),
)
