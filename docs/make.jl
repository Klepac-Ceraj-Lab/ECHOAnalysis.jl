using Documenter, ECHOAnalysis

makedocs(
    sitename="ECHO RESONANCE Microbiombe Paper 1",
    pages=[
    "Getting Started" => "index.md",
    "Dealing with Metadata" => "metadata_handling.md",
        ],
    format = Documenter.HTML(
            prettyurls = get(ENV, "CI", nothing) == "true"
    ),
)

deploydocs(
    repo="github.com/Klepac-Ceraj-Lab/echo_analysis"
)
